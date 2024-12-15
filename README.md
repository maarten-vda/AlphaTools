# AlphaTools
Tools for AlphaFold3

## Requirements

I've only tested it with one conda 4.14.0 environment running these versions

`python == 3.10.4`  
`pandas == 1.4.2`  
`plotly == 5.22.0`  
`matplotlib == 3.7.1`  
`numpy == 1.26.4`  
`biopython == 1.81`  
`requests == 2.29.0`  

A conda environment with these dependencies can be created with the following command.  
```conda create -n AlphaTools python=3.10.4 pandas=1.4.2 plotly=5.22.0 matplotlib=3.7.1 numpy=1.26.4 biopython=1.81 requests=2.29.0 -y

The script expects five models per fold folded with AlphaFold3. The input directory should be the unzipped output of an AlphaFold3 fold.  

# Running AlphaTools

## Interactive pAE Heatmap  

```python AlphaTools.py --input $INPUT_DIR --mode "pae"```   
This opens an interactive pAE heatmap in the default web browser, highlight areas to zoom in and double click to zoom out. This pAE heatmap is the mean pAE across the five folds, but this can be changed with the --combine flag.  
The --combine flag has `choices=['mean', 'range', 'max', 'min', 'median', 'std_dev']`, and the default value is 'mean'.  

To make a pAE heatmap with the minimum pAE across models, run the following command, and replace min with any other combine flag choices  
```python AlphaTools.py --input $INPUT_DIR --mode "pae" --combine "min"```  


## pLDDT plot

To make a pLDDT plot of all 5 models, run the command  
```python AlphaTools.py --input $INPUT_DIR --mode "plddt"```  

The --combine flag also works with pLDDT mode. For example the following command makes a pLDDT plot of the highest scores across models  
```python AlphaTools.py --input $INPUT_DIR --mode "plddt" --combine "max"```  


## Delta pAE  

This script can calculate the difference between the pAE heatmaps of two folds. This can be used to predict what the effect of missense mutations are on protein stability/binding.  
```python AlphaTools.py --input $INPUT_DIR --input2 $INPUT2_DIR --mode "delta_pae"```  

In this mode the --input2 is subtracted from --input1. `delta_pae` mode also supports the --combine flag.  
```python AlphaTools.py --input $INPUT_DIR --input2 $INPUT2_DIR --mode "delta_pae --combine "mean""```  
