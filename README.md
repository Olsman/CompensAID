CompensAID: An automated detection tool for spillover errors.


Overview
--------
CompensAID is an R package designed to automatically assesses which marker combinations exhibit signs of spillover errors in (spectral) flow cytometry data. 


Key Features
------------
- Calculate density-based cutoff-detection
- Perform automated gating
- Segment the positive population
- Calculate the secondary stain index (SSI) across each marker combination
- Output an SSI matrix and SSI informationd dataFrame. 


Installation
------------
devtools::install_github("Olsman/CompensAID")


Citation
------------
The CompensAID package integrates the SSI, as previously described:
Daniels, K. & Gardner, R. Secondary stain index. Memorial Sloan Kettering Cancer Center https://wi.mit.edu/sites/default/files/2021-05/20200504_Post-it_Secondary_Stain_Index_Final.pdf (2020).
