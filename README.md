# AlphaTools
Tools for AlphaFold3

## Requirements

I've only tested it with one environment running these versions

python == 3.10.4  
pandas == 1.4.2  
plotly == 5.22.0  
matplotlib == 3.7.1  
numpy == 1.26.4  
biopython == 1.81  
requests == 2.29.0  


## Running AlphaTools

# Interactive pAE Heatmap  

python AlphaTools.py --input $INPUT_DIR --mode "pae"  
This opens an interactive pAE heatmap in the default web browser, highlight areas to zoom in and double click to zoom out. This pAE heatmap is the mean pAE across the five folds, but this can be changed with the --combine flag.  
The --combine flag has choices=['mean', 'range', 'max', 'min', 'median', 'std_dev'], and the default value is 'mean'.  

To make a pAE heatmap with the minimum pAE, run the following command  
python AlphaTools.py --input $INPUT_DIR --mode "pae" --combine "min"


# pLDDT map

