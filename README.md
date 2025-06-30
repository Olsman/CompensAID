CompensAID: An Automated Detection Tool for Spillover Errors 


Overview
--------
CompensAID is an R package developed to automatically assess which marker combinations show signs of spillover/unmixing errors in both conventional and spectral flow cytometry data.


Key Features
--------
- Density-based cutoff detection
- Automated gating
- Positive population segmentation
- Calculation of the Secondary Stain Index (SSI) across all marker combinations
- Output of an SSI matrix and a detailed SSI information data frame


Installation
--------
Install the package using the following command in R:
devtools::install_github("Olsman/CompensAID")


Citation
--------
The CompensAID package incorporates the Secondary Stain Index (SSI) as described in:
Daniels, K., and Gardner, R. Secondary Stain Index. Memorial Sloan Kettering Cancer Center. https://wi.mit.edu/sites/default/files/2021-05/20200504_Post-it_Secondary_Stain_Index_Final.pdf (2020).
