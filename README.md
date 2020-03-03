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

## References
1. Kim, Y., Sidney, J., Pinilla, C. et al. Derivation of an amino acid similarity matrix for peptide:MHC binding and its application as a Bayesian prior. BMC Bioinformatics 10, 394 (2009). https://doi.org/10.1186/1471-2105-10-394
2. Robinson J, Maccari G, Marsh SGE, Walter L, Blokhuis J, Bimber B, Parham P, De Groot NG, Bontrop RE, Guethlein LA, and Hammond JA, KIR Nomenclature in non-human species, Immunogenetics (2018), in preparation
