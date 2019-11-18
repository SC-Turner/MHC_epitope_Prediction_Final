# MHC epitope Prediction
Screens MHC pocket sequences for functional homology at the pocket level

## Requirements
1. Ubuntu 16.04/18.04
2. R v3.6.1

## How to use
1. Create a new folder within /Inputfiles directory with your target MHC name
2. Place prob.txt (From I-Tasser/COACH) and target MHC .fasta into the new directory
3. Run the following code from the command line whilst in the /MHC_epitope_Prediction_Final directory
    Rscript --vanilla Extended_Pocket_Quantification.R <MHC/folder_name> <Pocket_Letter_(B/F)>
4. Output .Rdata and .png results will export to /Outputfiles Directory

*Example input files used for the current project have been included in /Inputfiles
