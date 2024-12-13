import argparse
import json
import pandas as pd
from collections import Counter
import plotly.graph_objects as go
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
import os
import glob
from Bio.Align import PairwiseAligner
import requests 


'''This section of the script is for helper functions to be used in the main scripts'''

def three_letter_to_one_letter(residue):
    # This function converts a three-letter amino acid code to a one-letter code
    residue_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return residue_dict[residue]

def cif_to_chains_sequence_dict(cif_filepath):
    ##Get the cif file from the corresponding json file
    chains_sequence_dict = {}
    with open(cif_filepath, 'r') as file:
        # Initialize variables
        prev_res_id = 0
        chain_id = None
        res_id = -1
        res = None
        for line in file:
            if line.startswith("ATOM"):
                line = line.split()
                atom_id = line[1]
                res_id = int(line[8])  # Convert to integer
                res = line[5]
                chain_id = line[6]
                if chain_id not in chains_sequence_dict:
                    chains_sequence_dict[chain_id] = ""
                if res_id != prev_res_id:
                    prev_res_id = res_id
                    chains_sequence_dict[chain_id] += three_letter_to_one_letter(res)
    return chains_sequence_dict

def cif_to_residue_atoms_dict(filepath):
    ##Get the cif file from the corresponding json file
    filename_prefix = filepath.rsplit('_', 3)[0]
    model_number = filepath.rsplit('_', 1)[1].split('.')[0]
    cif_filepath = f"{filename_prefix}_model_{model_number}.cif"
    residue_atoms_dict = {}
    with open(cif_filepath, 'r') as file:
        # Initialize variables
        prev_chain_id = None
        prev_chain_max_res_id = 0
        prev_res_id = 0

        for line in file:
            if line.startswith("ATOM"):
                line = line.split()
                atom_id = line[1]
                res_id = int(line[8])  # Convert to integer
                chain_id = line[6]

                # Update prev_chain_max_res_id if chain_id changes
                if chain_id != prev_chain_id:
                    prev_chain_id = chain_id
                    prev_chain_max_res_id = prev_res_id

                # Calculate adjusted_res_id
                adjusted_res_id = res_id + prev_chain_max_res_id

                # Update residue_atoms_dict
                if adjusted_res_id not in residue_atoms_dict:
                    residue_atoms_dict[adjusted_res_id] = 1
                else:
                    residue_atoms_dict[adjusted_res_id] += 1

                # Update prev_res_id
                prev_res_id = res_id + prev_chain_max_res_id

    return residue_atoms_dict

def pair_inputs(fold_selection_string_pairs):
    # Initialize an empty list to store the pairs
    # This function is needed to have an unlimited number of inputs and different formats for each input
    pairs = []

    # Iterate over the input list in steps of 2
    for i in range(0, len(fold_selection_string_pairs), 2):
        # Append the current pair to the list
        pairs.append((fold_selection_string_pairs[i], fold_selection_string_pairs[i + 1]))

    return pairs

def parse_ranges(input_list):
    # This function uses regex to convert the chain selection string to something machine readable
    # Initialize an empty dictionary to hold results
    result = {}

    # Regular expression to match entries with ranges
    range_pattern = re.compile(r"([A-Z])(\d+)-(\d+)")
    # Regular expression to match single letter entries without ranges
    letter_pattern = re.compile(r"^[A-Z]$")

    for entry in input_list:
        # Check if entry matches the range pattern
        if match := range_pattern.match(entry):
            letter, start, end = match.groups()
            # Add the range to the corresponding letter key
            if letter not in result:
                result[letter] = []
            result[letter].append((int(start), int(end)))
        # Check if entry is a single letter without ranges
        elif letter_pattern.match(entry):
            letter = entry
            if letter not in result:
                result[letter] = []

    return result

def translate_selection_string(selection_string):
    # This function converts the whole selection string to machine readable format
    selection_dict = {}
    split_string = selection_string.split('/')
    selection_dict['residues'] = parse_ranges(split_string[0].split(';'))
    selection_dict['models'] = list(split_string[1])
    selection_dict['method'] = str(split_string[2])
    return(selection_dict)

def get_chain_boundaries(input_file):
    with open(input_file, 'r') as file:
        data = json.load(file)
        token_chain_ids = data.get("token_chain_ids", [])
        element_counts = Counter(token_chain_ids)
        return element_counts

def load_json_plddt_data(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
        plddt_scores = data.get("atom_plddts", [])
        # Function that gets which atoms correspond to which residues
        residue_atoms_dict = cif_to_residue_atoms_dict(filepath)
        per_residue_plddt_scores = []
        prev_atom_count = 0
        # Divide the summed PLDDT scores by the number of atoms in each residue
        for residue_id, atom_count in residue_atoms_dict.items():
            residue_plddt_score = sum(plddt_scores[prev_atom_count:prev_atom_count + atom_count]) / atom_count
            per_residue_plddt_scores.append(residue_plddt_score)
            prev_atom_count += atom_count
        return per_residue_plddt_scores

def align_folds(paired_inputs):
    '''This part extracts the sequence of each chain from the corresponding CIF files. 
    It then aligns the folds with Bio.Align.PairwiseAligner and puts the scores in a DF, 
    offset_df which shows which sequences are the same binarily (is vs None) as well as 
    the relative offset position of fold x to fold y. If there are repeated proteins, 
    there is logic to pair the chains alphabetically for chains in the selection string'''
    chain_sequence_dict = {}
    for i in paired_inputs:
        directory = i[0]
        # Use glob to find the file matching the pattern
        cif_filepath = [os.path.join(directory, filename) 
                  for filename in os.listdir(directory) 
                  if filename.endswith("model_0.cif")][0]

        chain_sequence_dict[directory] = cif_to_chains_sequence_dict(cif_filepath)
        
    ### Initialize an empty dataframe with the dimensions of the number of chains
    alignment_df = pd.DataFrame(np.empty((len(list(chain_sequence_dict.values())[0]), len(list(chain_sequence_dict.values())[1])), dtype=object))  ## This df keeps the alignment score
    relative_pos_df = pd.DataFrame(np.empty((len(list(chain_sequence_dict.values())[0]), len(list(chain_sequence_dict.values())[1])), dtype=object))  ## This df keeps the relative positions of the aligned sequences (raw PairwiseAligner output)
    offset_df = pd.DataFrame(np.empty((len(list(chain_sequence_dict.values())[0]), len(list(chain_sequence_dict.values())[1])), dtype=object)) ## This df keeps the start and end positions of the seq2 relative to seq1

    fold1 = list(chain_sequence_dict.keys())[0]
    fold2 = list(chain_sequence_dict.keys())[1]

    ## Assign the header row and column as the chain IDs
    alignment_df.index = chain_sequence_dict[fold1].keys()  
    alignment_df.columns = chain_sequence_dict[fold2].keys() 
    relative_pos_df.index = chain_sequence_dict[fold1].keys()
    relative_pos_df.columns = chain_sequence_dict[fold2].keys()
    offset_df.index = chain_sequence_dict[fold1].keys()
    offset_df.columns = chain_sequence_dict[fold2].keys()


    for chain1, seq1 in (chain_sequence_dict[fold1]).items():
        for chain2, seq2 in (chain_sequence_dict[fold2]).items():
            relative_pos_df.loc[chain1][chain2], alignment_df.loc[chain1][chain2] = align_sequences_alt_method(seq1, seq2) ## Align the sequences and save the score in the dataframe

    alignment_df = alignment_df.astype(float) ## Convert the dataframe to float
    alignment_df[alignment_df < 0.90] = 0 ## Set the threshold for alignment to 0.95, to account for mutations and insertions where the sequence wont match 100%
    alignment_df[alignment_df >= 0.90] = 1 ## Binarise the matrix
    relative_pos_df[alignment_df == 0] = None ## Get rid of the relative positions of unaligned sequences



    ## Logic for getting the offset positions
    ### THIS SHOULD WORK BUT I HAVENT TESTED IT YET ON ANY FOLDS JUST RANDOM SEQUENCES
    for index, row in relative_pos_df.iterrows():
        for column in row.index:
            if relative_pos_df.loc[index][column] is not None:
                alignment = relative_pos_df.loc[index][column]
                start1 = alignment[0][0][0]
                start2 = alignment[1][0][0]
                offset_df.loc[index][column] = start1 - start2
            else:
                offset_df.loc[index][column] = None 


    '''This section is the logic for choosing which repeated proteins are overlapped in the pLDDT plot. 
    This depends on what chains are chosen to be shown, and otherwise aligns the chains in alphbetical order'''

    ### This is the logic for identifying repeat proteins, makes repeat_protein_dict
        ### First identifies sites with more than 1 1 in the row or column
        ### Takes column and row IDs where theres more than 1 1 and adds to a dictionary

    fold1_total = []
    fold2_total = []
    repeat_protein_dict = {} 

    for row_index, row in alignment_df.iterrows():
        row_one_count = (row == 1).sum()  # Count occurrences of 1 in the row
        for column in alignment_df.columns:
            column_one_count = (offset_df[column] == 1).sum()  # Count occurrences of 1 in the column
            if column not in fold2_total and row_index not in fold1_total:
                if column_one_count > 1 or row_one_count > 1:
                    row_indices = alignment_df.index[alignment_df.loc[row_index] == 1].tolist()
                    column_indices = alignment_df.columns[alignment_df[row_index] == 1].tolist()
                    repeat_protein_dict[row_index] = {f"{fold1}": row_indices, f"{fold2}": column_indices}
                    fold1_total.extend(row_indices)
                    fold2_total.extend(column_indices)

    ## This part figures out which chains are shown in which folds
    which_chains = {}
    for i in paired_inputs:
        filepath = i[0]
        which_chains[filepath] = []
        for key, value in (i[1]['residues']).items():
            which_chains[filepath].append(key)       

    # This removes everything in repeat_protein_dict not shown in the selection string
    for repeat_name, repeats in repeat_protein_dict.items():
        for key, value in repeats.items():
            to_remove = []  # Temporary list to store items to remove
            for i in value:
                if i not in which_chains[key]:
                    to_remove.append(i)  # Mark `i` for removal
            # Remove items after the inner loop
            for i in to_remove:
                repeat_protein_dict[repeat_name][key].remove(i)

    ### This part iterates through repeat_protein_dict, then iterates whichever list is shorter (alphabetically since its already sorted that way)
    ### Then it sets everything in the same row and column except the alphabetically corresponding cell

    for repeat_name, repeats in repeat_protein_dict.items():
        shortest_len = min(len(value) for value in repeats.values())
        paired_repeats = {key: sorted(value)[:shortest_len] for key, value in repeats.items()}
        chains1 = list(repeats.values())[0]
        chains2 = list(repeats.values())[1]
        ## Here need to make chains1 and chains2 so I can get the nth term
        for fold, chains in paired_repeats.items():
            if fold == fold1:
                ## Here I need to iterate n through the length of chains
                ## Set everything in the offset_df column with the same label as pos n in chains to none 
                # except the row with the label in pos n of fold2
                for i in range(len(chains1)):
                    fold1_chain = chains1[i]
                    fold2_chain = chains2[i]
                    # Set all values in the specified column to None except for the intersecting value
                    offset_df[fold1_chain] = offset_df[fold1_chain].where(offset_df.index == fold2_chain, None)

                    # Set all values in the specified row to None except for the intersecting value
                    offset_df.loc[fold2_chain] = offset_df.loc[fold2_chain].where(offset_df.columns == fold1_chain, None)

            ## Row = fold1, column = fold2

    return (offset_df)

    # Perform global alignment
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    return(alignments)

def align_sequences_alt_method(sequence1, sequence2):
    # Perform global alignment
    aligner = PairwiseAligner()
    aligner.gap_score = 0  ### I set the gap score to 0 to only account for mismatches, if I am trying to align a truncated vs full-length protein
    aligner.mismatch_score = -1
    aligner.match_score = 1
    alignments = aligner.align(sequence1, sequence2)
    return(alignments[0].aligned, alignments[0].score/min(len(sequence1), len(sequence2))) ### Score divided by the min length between both sequences so its normalised, max alignment score is same as min length

def select_files(paired_inputs):
    # This function chooses which files to load based on the selection string
    # Input is paired_inputs string with translated selection strings, output is dictionary with directory as key and list of files as value
    filepaths_dict = {}
    for i in paired_inputs:
        dir = i[0]
        models = i[1]['models']
        files = []
        for model in models:
            # Create the pattern for each model
            pattern = os.path.join(dir, f"*full_data_{model}.json")
            
            # Use glob to search for matching files
            matching_files = glob.glob(pattern)
            
            # Extend the list with found files
            files.extend(matching_files)
        filepaths_dict[dir] = files
    return filepaths_dict

def convert_offset_selection_string_to_absolute(paired_inputs, boundaries_dict):
    ### This function will take the offset_df and the boundaries_dict and convert the offset to the absolute residue ID

    offset_df = align_folds(paired_inputs)

    fold1 = paired_inputs[0][0]
    fold2 = paired_inputs[1][0]

    ## Iterate through boundaries dict and add them up to make the absolute pos dict
    absolute_pos_dict = {}
    for fold, boundaries in boundaries_dict.items():
        total = 0
        absolute_pos_dict[fold] = {}
        for i, j in boundaries.items():
            total += j
            absolute_pos_dict[fold][i] = total

    ## Relate the absolute position in pLDDT scores to the offset
    pos_offset_dict = {}
    for fold, boundaries in absolute_pos_dict.items():
        prev_pos = 0
        cumulative_offset = 0
        pos_offset_dict[fold] = {}
        for chain, absolute_pos in boundaries.items():
            if fold == fold1:
                ## Dont offset the first fold since we are oriented to the first fold
                pos_offset_dict[fold][chain] = {"plddt_start": prev_pos + 1, "plddt_end": absolute_pos, "offset": 0}
            elif fold == fold2:
                offset = offset_df.loc[chain].dropna().iloc[0]
                cumulative_offset += offset
                pos_offset_dict[fold][chain] = {"plddt_start": prev_pos + 1, "plddt_end": absolute_pos, "offset": cumulative_offset}
            prev_pos = absolute_pos

    ### Add processing here so that if theres a negative offset it offsets the first one

    ### This part combines which residues are shown in the selection string, and the offset df
    ### to make a dictionary relating positions in each chain between the json data and plot
    ### Dictionary with fold: {chain: {plddt_start: absolute_start, plddt_end: absolute_end, plot_start: plot_start, plot_end: plot_end}}

    json_plot_association = {}

    for fold, value in pos_offset_dict.items():
        for t in paired_inputs:
            if t[0] == fold:
                selection_string = t[1]
                residues = selection_string['residues']
        json_plot_association[fold] = {}
        prev_json_end = 0
        for chain, pos_dict in value.items():

            json_start = pos_dict['plddt_start'] - 1 ## -1 because the selection string is 1 oriented and python is 0
            json_end = pos_dict['plddt_end'] - 1

            offset = pos_dict['offset'] ## This is the absolute offset at that residue so is cumulative
            try:
                window_list = residues[chain]
                if len(window_list) == 0: # If window (which is the per-chain residues to show in selection string) has length of 0, means show everything
                    window_json_start = json_start
                    window_json_end = json_end
                elif len(window_list) > 0:
                    window = window_list[0] ##### I have set this to the first element in the list, assuming it is a list of 1. Need to add more processing to handle aligning length changes in the middle of a protein
                    window_start = window[0] - 1 ## -1 because the selection string is 1 oriented and python is 0
                    window_end = window[1] - 1 
                    window_len = window_end - window_start
                    window_json_start = json_start - window_start ## Window_json_start means the start pos of the chain in the json corrected for the selection string 
                    window_json_end = window_json_start + window_len
                    if fold == fold2:
                        ()
                    elif fold == fold1:
                        ()
                        #plot_start = 

            ### How do I adjust here plot_start and plot_end for the selection string, offset, and other fold
            ### Can always have relative pos between folds have negative start and gets adjusted
            ### How do I preserve the alignment and still trim?

            except KeyError: # If a chain is not found in the selection string do not put it in json plot association
                continue 
    return()


def combine_datasets(filepaths, combination_method='mean'):
    """
    Combines datasets from multiple JSON files using the specified combination method.

    Parameters:
    - filepaths (list of str): List of file paths to the JSON files.
    - combination_method (str): Method to combine the data. Can be one of ['mean', 'range', 'max', 'min', 'median', 'std_dev'].

    Returns:
    - pd.DataFrame: The combined DataFrame based on the specified combination method.
    """
    
    # Initialize a list to store DataFrames
    data_frames = []

    # Load JSON data from each file
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            data = json.load(file)
            predicted_aligned_error = data.get("pae", [])
            df = pd.DataFrame(predicted_aligned_error)
            data_frames.append(df)

    # Check if there are any data frames to combine
    if not data_frames:
        raise ValueError("No data to combine. Please check the filepaths and the data format.")

    # Combine the DataFrames based on the specified method
    if combination_method == 'mean':
        combined_df = pd.concat(data_frames, axis=0).groupby(level=0).mean()  # Calculate the mean of all rows
    elif combination_method == 'range':
        combined_df = pd.concat(data_frames, axis=0).groupby(level=0).apply(lambda x: x.max() - x.min())  # Calculate the range (max - min)
    elif combination_method == 'max':
        combined_df = pd.concat(data_frames, axis=0).groupby(level=0).max()  # Take the max value
    elif combination_method == 'min':
        combined_df = pd.concat(data_frames, axis=0).groupby(level=0).min()  # Take the min value
    elif combination_method == 'median':
        combined_df = pd.concat(data_frames, axis=0).groupby(level=0).median()  # Calculate the median
    elif combination_method == 'std_dev':
        combined_df = pd.concat(data_frames, axis=0).groupby(level=0).std()  # Calculate the standard deviation
    else:
        raise ValueError(f"Unsupported combination method: {combination_method}")

    return combined_df

def get_contacts(paired_inputs):
    return()

# Function to extract unique numbers and calculate their mean
def extract_and_average(confidence_value):
    # Extract numbers from the string
    numbers = [float(num) for num in re.findall(r'\d+\.\d+|\d+', confidence_value)]

    # Return the mean of unique numbers if there are any, otherwise return the original value
    if numbers:
        return np.mean(np.unique(numbers))
    return confidence_value  # In case no numbers are found, keep the original value


def get_complex_data(uniprot_id):

    intact_filepath = "~/Desktop/human_complex_portal.tsv"
    # Specify the columns to load
    columns_to_load = ["#Complex ac", 'Recommended name', "Identifiers (and stoichiometry) of molecules in complex", 'Expanded participant list']

    df = pd.read_csv(intact_filepath, sep='\t', usecols=columns_to_load)

    df.rename(columns={"Identifiers (and stoichiometry) of molecules in complex": 'Participant list'}, inplace=True)

    filtered_df = df[df['Expanded participant list'].str.contains(uniprot_id, case=False, na=False)]

    rows_to_drop = []

    # This part converts the participants list string to a dictionary, and drops instances where the number of times a protein appears is 0 (RNA stuff)
    for index, row in filtered_df.iterrows():
        split_list = row['Expanded participant list'].split('|')
        participant_dict = {}
        skip_row = False
        for i in split_list:
            value_inside_brackets = i.split('(')[1].split(')')[0]
            participant_id = i.split('(')[0]
            if value_inside_brackets != str(0):
                participant_dict[participant_id] = value_inside_brackets
            if value_inside_brackets == str(0) and participant_id == uniprot_id:
                skip_row = True
                rows_to_drop.append(index)
                break
        if skip_row:
            continue
        filtered_df.at[index, 'Expanded participant list'] = participant_dict

    filtered_df = filtered_df.drop(rows_to_drop)

    print(filtered_df)

    return()

'''This section of the script is for the data extraction functionalities of AlphaTools'''
def get_value_from_file(input_file, x, y):
    with open(input_file, 'r') as file:
        data = json.load(file)
        predicted_aligned_error = data.get("pae", [])
        if 1 <= x <= len(predicted_aligned_error):
            if 1 <= y <= len(predicted_aligned_error):
                return predicted_aligned_error[x-1][y-1]
            else:
                print(f"Element y ({y}) is out of range in list number {x}")
        else:
            print(f"Element x ({x}) is out of range")

def extract_metrics(input_dir):
    '''This script extracts several quality metrics from the JSON files in the input directory
    It saves three TSV files, one detailing all the contacts and one detailling the interfaces
    The contacts file has rows corresponding to each unique contact. The contacts file has these columns:

    - chain1: The ID of the first chain in the contact
    - chain2: The ID of the second chain in the contact
    - residue1: The residue number in the first chain
    - residue2: The residue number in the second chain
    - res1_code: The one-letter code of the residue in the first chain
    - res2_code: The one-letter code of the residue in the second chain
    - min_dist: The minimum distance between the two residues
    - mean_dist: The mean distance between the two residues
    - dist_std_dev: The standard deviation of distances between the two residues
    - num_atom_pairs: The number of atom pairs between the two residues that classify as a contact
    - num_h_bonds: The number of hydrogen bonds between the two residues
    - num_salt_bridges: The number of salt bridges between the two residues
    - num_pi_stacks: The number of pi-stacking interactions between the two residues
    - model_count: The number of models in which the contact is present
    - models: A list of models in which the contact is present
    - min_pae: The minimum predicted aligned error between the two residues
    - mean_pae: The mean predicted aligned error between the two residues
    - pae_std_dev: The standard deviation of predicted aligned errors between the two residues
    - min_plddt: The minimum pLDDT score between the two residues
    - mean_plddt: The mean pLDDT score between the two residues
    - plddt_std_dev: The standard deviation of pLDDT scores between the two residues
    - min_contact_prob: The minimum contact probability between the two residues
    - mean_contact_prob: The mean contact probability between the two residues
    - contact_prob_std_dev: The standard deviation of contact probabilities between the two residues

        I could potentially add something with ∆pLDDT to the fold by itself

    The interfaces file has rows corresponding to each unique interface. The interfaces file has these columns:
    - interface: which chains, formatted as {chain1}:{chain2}
    - averages for everything in the contacts file
    - num_contacts: The number of unique residue contacts in the interface
    - num_contacts_in_all_models: The number of unique residue contacts present in all models
    - num_contacts_min_4_models: The number of unique residue contacts present in at least 4 models
    - num_contacts_min_3_models: The number of unique residue contacts present in at least 3 models
    - num_contacts_min_2_models: The number of unique residue contacts present in at least 2 models
    - chain_pair_iptm: ipTM score for the chain pair (in summary_confidences_{n}.json)
    - chain_pair_pae_min: Minimum pAE score for the chain pair
    
    I also want to include metrics relating to SASA, buried surface area, interaction energy

    Then I need another TSV, with information on the complex as a whole. 
    This should have information about what chains are complex 1 and what chains are complex 2.
    Then it can combine the metrics from the interfaces file to give a summary of the predicted interaction as a whole.
    Since this file is per complex, it can include other metrics from summary_confidences_{n}.json like pTM, fraction disordered, has_clash, etc.
    In this file I can put additional features for ML, like the protein feature libraries

    '''
    
    
    return()

'''This section of the script is for the graphing functionalities of AlphaTools'''

def delta_pae(input_dir, input_dir2, combination_method):
    # Initialize an empty list to store the file paths
    filepaths = []
    filepaths2 = []

    # Loop through all files in the directory
    for filename in os.listdir(input_dir):
        # Use regular expression to match the pattern
        if re.match(r'.*full_data_\d+\.json$', filename):
            # If the filename matches the pattern, add it to the list
            filepaths.append(os.path.join(input_dir, filename))

    for filename in os.listdir(input_dir2):
        if re.match(r'.*full_data_\d+\.json$', filename):
            filepaths2.append(os.path.join(input_dir2, filename))

    # Combine the datasets using the specified method
    df1 = combine_datasets(filepaths, combination_method)
    df2 = combine_datasets(filepaths2, combination_method)

    df = df1 - df2

    # Adjust the index and columns to start from 1
    df.index = df.index + 1
    df.columns = df.columns + 1
    
    # Get cumulative boundaries
    boundaries = get_chain_boundaries(filepaths[0])
    name1 = filepaths[0].rsplit('/', 2)[-2].split('full_data')[0].strip('/')

    # Calculate cumulative boundaries
    cumulative_boundaries = {}
    lower_bound = 1
    for key in boundaries:
        upper_bound = lower_bound + boundaries[key] - 1
        cumulative_boundaries[key] = (lower_bound, upper_bound)
        lower_bound = upper_bound + 1

    #colorscale that matches alphafold-multimer
    #red is #ff0808
    #blue is #0404ff

    # Custom colorscale: red (#ff0808) at 0, white at 0.5, and blue (#0404ff) at 1
    af_rdbu = [
        [0, 'rgb(4, 4, 255)'],   # Blue
        [0.5, 'rgb(255, 255, 255)'],  # White
        [1, 'rgb(255, 8, 8)']  # Red

    ]

    # Create an interactive heatmap using Plotly with WebGL
    fig = go.Figure(data=go.Heatmap(
        z=df.values,
        x=df.columns,
        y=df.index,
        colorscale=af_rdbu,
        colorbar=dict(title="PAE"),
        zsmooth='best',  # Smooth zooming
    ))

    # Add boundary lines
    for key, (start, end) in cumulative_boundaries.items():
        fig.add_shape(type="line",
                      x0=start-0.5, y0=0.5, x1=start-0.5, y1=len(df)+0.5,
                      line=dict(color="black", width=2))
        fig.add_shape(type="line",
                      x0=0.5, y0=start-0.5, x1=len(df)+0.5, y1=start-0.5,
                      line=dict(color="black", width=2))

    fig.update_layout(
        title="Predicted Aligned Error (PAE) Heatmap for " + name1 + f" ({combination_method} of 5 models)",
        xaxis=dict(
            title="Residue 2",
            tickmode="array",  # Set tick mode to "array"
            tickvals=[(cumulative_boundaries[key][0] + cumulative_boundaries[key][1]) / 2 for key in list(boundaries.keys())],
            ticktext=list(boundaries.keys()),
        ),
        yaxis=dict(
            title="Residue 1 (error aligned to this residue)",
            tickmode="array",  # Set tick mode to "array"
            tickvals=[(cumulative_boundaries[key][0] + cumulative_boundaries[key][1]) / 2 for key in list(boundaries.keys())],
            ticktext=list(boundaries.keys()),
            autorange="reversed"  # Ensure the y-axis starts from the top
        )
    
    )

    # Show the heatmap
    fig.show()
    pass

def avg_heatmap(input_dir, combination_method):
    
    # Initialize an empty list to store the file paths
    filepaths = []

    # Loop through all files in the directory
    for filename in os.listdir(input_dir):
        # Use regular expression to match the pattern
        if re.match(r'.*full_data_\d+\.json$', filename):
            # If the filename matches the pattern, add it to the list
            filepaths.append(os.path.join(input_dir, filename))

    # Combine the datasets using the specified method
    df = combine_datasets(filepaths, combination_method)
    
    # Adjust the index and columns to start from 1
    df.index = df.index + 1
    df.columns = df.columns + 1
    
    # Get cumulative boundaries
    boundaries = get_chain_boundaries(filepaths[0])
    name1 = filepaths[0].rsplit('/', 2)[-2].split('full_data')[0].strip('/')

    # Calculate cumulative boundaries
    cumulative_boundaries = {}
    lower_bound = 1
    for key in boundaries:
        upper_bound = lower_bound + boundaries[key] - 1
        cumulative_boundaries[key] = (lower_bound, upper_bound)
        lower_bound = upper_bound + 1

    #colorscale that matches alphafold-multimer
    #red is #ff0808
    #blue is #0404ff

    # Custom colorscale: red (#ff0808) at 0, white at 0.5, and blue (#0404ff) at 1
    af_rdbu = [
        [0, 'rgb(4, 4, 255)'],   # Blue
        [0.5, 'rgb(255, 255, 255)'],  # White
        [1, 'rgb(255, 8, 8)']  # Red

    ]

    # Create an interactive heatmap using Plotly with WebGL
    fig = go.Figure(data=go.Heatmap(
        z=df.values,
        x=df.columns,
        y=df.index,
        colorscale=af_rdbu,
        colorbar=dict(title="PAE"),
        zsmooth='best',  # Smooth zooming
    ))

    # Add boundary lines
    for key, (start, end) in cumulative_boundaries.items():
        fig.add_shape(type="line",
                      x0=start-0.5, y0=0.5, x1=start-0.5, y1=len(df)+0.5,
                      line=dict(color="black", width=2))
        fig.add_shape(type="line",
                      x0=0.5, y0=start-0.5, x1=len(df)+0.5, y1=start-0.5,
                      line=dict(color="black", width=2))

    fig.update_layout(
        title="Predicted Aligned Error (PAE) Heatmap for " + name1 + f" ({combination_method} of 5 models)",
        xaxis=dict(
            title="Residue 2",
            tickmode="array",  # Set tick mode to "array"
            tickvals=[(cumulative_boundaries[key][0] + cumulative_boundaries[key][1]) / 2 for key in list(boundaries.keys())],
            ticktext=list(boundaries.keys()),
        ),
        yaxis=dict(
            title="Residue 1 (error aligned to this residue)",
            tickmode="array",  # Set tick mode to "array"
            tickvals=[(cumulative_boundaries[key][0] + cumulative_boundaries[key][1]) / 2 for key in list(boundaries.keys())],
            ticktext=list(boundaries.keys()),
            autorange="reversed"  # Ensure the y-axis starts from the top
        )
    
    )

    # Show the heatmap
    fig.show()
    pass

def plddt_per_residue(input_dir, combination_method):
    '''Write a function that takes a json file and returns a PLDDT plot.'''
    # Initialize an empty list to store the file paths
    filepaths = []

    # Loop through all files in the directory
    for filename in os.listdir(input_dir):
        # Use regular expression to match the pattern
        if re.match(r'.*full_data_\d+\.json$', filename):
            # If the filename matches the pattern, add it to the list
            filepaths.append(os.path.join(input_dir, filename))

    # Initialize a list to store all PLDDT scores
    all_plddt_scores = []

    # Load JSON data from each file
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            data = json.load(file)
            plddt_scores = data.get("atom_plddts", [])
            #Function that gets which atoms correspond to which residues
            residue_atoms_dict = cif_to_residue_atoms_dict(filepath)
            per_residue_plddt_scores = []
            prev_atom_count = 0
            #Divide the summed PLDDT scores by the number of atoms in each residue
            for residue_id, atom_count in residue_atoms_dict.items():
                if combination_method == 'mean':
                    residue_plddt_score = sum(plddt_scores[prev_atom_count:prev_atom_count + atom_count]) / atom_count
                elif combination_method == 'max':
                    residue_plddt_score = max(plddt_scores[prev_atom_count:prev_atom_count + atom_count])
                elif combination_method == 'min':
                    residue_plddt_score = min(plddt_scores[prev_atom_count:prev_atom_count + atom_count])
                elif combination_method == 'median':
                    residue_plddt_score = np.median(plddt_scores[prev_atom_count:prev_atom_count + atom_count])
                elif combination_method == 'std_dev':
                    residue_plddt_score = np.std(plddt_scores[prev_atom_count:prev_atom_count + atom_count])
                elif combination_method == 'range':
                    residue_plddt_score = max(plddt_scores[prev_atom_count:prev_atom_count + atom_count]) - min(plddt_scores[prev_atom_count:prev_atom_count + atom_count])
                per_residue_plddt_scores.append(residue_plddt_score)
                prev_atom_count += atom_count
            all_plddt_scores.append(per_residue_plddt_scores)

    # Get cumulative boundaries
    boundaries = get_chain_boundaries(filepaths[0])

    # Calculate cumulative boundaries
    cumulative_boundaries = {}
    lower_bound = 1
    for key in boundaries:
        upper_bound = lower_bound + boundaries[key] - 1
        cumulative_boundaries[key] = (lower_bound, upper_bound)
        lower_bound = upper_bound + 1

    

    # Plot the PLDDT scores
    plt.figure(figsize=(10, 6))
    
        
    # Add boundary lines
    for key, (start, end) in cumulative_boundaries.items():
        # Add vertical lines
        plt.axvline(x=start - 0.5, color='black', linewidth=1.5)
        
        # Add horizontal lines
        plt.axhline(y=start - 0.5, color='black', linewidth=1.5)

    min_x = float('inf')
    max_x = float('-inf')
    min_y = float('inf')
    max_y = float('-inf')

    for i, plddt_scores in enumerate(all_plddt_scores):
        x_values = range(1, len(plddt_scores) + 1)
        plt.plot(x_values, plddt_scores, label=f'Model {i + 1}', linewidth=0.5, alpha=0.7)
                # Update min and max values
        min_x = min(min_x, min(x_values))
        max_x = max(max_x, max(x_values))
        min_y = min(min_y, min(plddt_scores))
        max_y = max(max_y, max(plddt_scores))
    

    plt.xlabel('Residue Number')
    plt.ylabel('PLDDT Score')
    plt.title('PLDDT Scores Across 5 Models')

    # Set the axis limits
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)

    plt.legend()
    plt.grid(False)
    plt.show()

def overlap_plddts(paired_inputs):
    filepaths_dict = select_files(paired_inputs)
    all_plddt_scores_dict = {}
    boundaries_dict = {}
    ###Iterate through the filepaths dictionary and load the json files into plddt scores dictionary
    for key, value in filepaths_dict.items():
        all_plddt_scores_dict[key] = []
        #Load boundaries once for each fold since its the same across folds, boundaries needed to find start and end of each chain for alignment
        boundaries_dict[key] = get_chain_boundaries(value[0])
        for i in value:
            ## Add functionality here to get absolute ID of resis to include
            
            all_plddt_scores_dict[key].append(load_json_plddt_data(i))  
            

    ## Here I need to align the different folds and find which residues to plot and where   
    convert_offset_selection_string_to_absolute(paired_inputs, boundaries_dict)

    ## Processing to combine pLDDT scores based on selection string

    ## Processing to select + plot different residues


    #### I dont need to make pLDDT per residue script, mono-pLDDT should work in this function and I can adjust names later
    ### Actually, does mono-pLDDT work in this function? Somehow I think trying to align just one fold will fail since theres never any self-comparison
    ### So I should figure out how to plot the pLDDT of each chain individually, offsetting by the offset and trimming plotted bit by selection string
    ### I think I can just copy the function for mono-pLDDT and remove the align_folds part, much easier to write one function for n=1 and one for n>1 rather than a very complicated function for all
    ### And I think I can achieve this by literally plotting out start-end for each chain, might have to modify above to include protein length


    return()


'''This section of the script is for utility functions of AlphaTools'''
def fasta_to_af3_json_input(fasta, output):
    ()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve value from a dictionary in a JSON file")
    parser.add_argument('fold_selection_string_pairs', metavar='fold_selection', type=str, nargs='*', help="Pairs of inputs and format strings, e.g., input1 'input1_format' input2 'input2_format'.")
    parser.add_argument('--mode', type=str, help="The mode to run the script in", default="none", required=True)
    # Add mutually exclusive flag for backend mode
    parser.add_argument('--backend', action='store_true', help="Enable backend mode, requires --output argument")
    # Add output argument without requiring it initially
    parser.add_argument('--output', type=str, help="Filepath for output png, required if --backend is set")
    parser.add_argument('--input', type=str, help="Input directory filepath or JSON file", required=False)
    parser.add_argument('--input2', type=str, help="Second input directory filepath or JSON file", required=False)
    parser.add_argument('--x', type=int, help="X-coordinate to retrieve value from", required=False)
    parser.add_argument('--y', type=int, help="Y-coordinate to retrieve value from", required=False)
    parser.add_argument('--combine', type=str, help="Way to combine two files", required=False, default='mean', choices=['mean', 'range', 'max', 'min', 'median', 'std_dev'])
    parser.add_argument('--fasta', type=str, help="Fasta file to convert to AlphaFold3 input JSON format", required=False)

    # Parse arguments and translate selection strings
    args = parser.parse_args()
    if args.fold_selection_string_pairs:
        paired_inputs = pair_inputs(args.fold_selection_string_pairs)
        for i in range(len(paired_inputs)):
            paired_inputs[i] = (paired_inputs[i][0], translate_selection_string(paired_inputs[i][1]))

    # Conditional logic to enforce --output only if --backend is selected
    if args.backend and not args.output:
        parser.error("--output is required when --backend is specified.")
    elif not args.backend and args.output:
        parser.error("--output can only be used with --backend mode.")

    print(get_complex_data("Q5TBB1"))

    if args.mode == "pae_value":
        value = get_value_from_file(args.input, args.x, args.y)
        print(f"The value at element {args.x},{args.y} is: {value}")
    if args.mode == "pae":
        avg_heatmap(args.input, args.combine)
    if args.mode == "delta_pae":
        delta_pae(args.input, args.input2, args.combine)
    if args.mode == "plddt":
        plddt_per_residue(args.input, args.combine)
    if args.mode == "overlap_plddts":
        overlap_plddts(paired_inputs)
    if args.mode == "make_input_json":
        fasta_to_af3_json_input(args.fasta, args.output)
    if args.mode == "get_interactors":
        get_complex_data(args.input)