# eht-imaging Sgr A* Stokes I Static Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** May 12, 2022

**Primary Reference:** The Event Horizon Telescope Collaboration, et al. 2022, ApJL (Sgr A* Paper III)

**Data Product Code:** [2022-D01-02](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

The pipeline reconstructs an image from uvfits files simultaneously released in the EHT website (data release ID: 2022-D01-01) using the python-interfaced imaging package eht-imaging (Chael et al. 2018).

To run the pipeline, specify the input uvfits data file, the file containing the imaging parameters, and the combination to be used. Multiple bands can be specified separately. Additional flags control the output, which is only the reconstructed image as a FITS file by default.

Example call:

    python eht-imaging_pipeline.py -i ../../data/uvfits_norm/ER6_SGRA_2017_097_lo_hops_netcal-LMTcal-norm_StokesI.uvfits -i2 ../../data/uvfits_norm/ER6_SGRA_2017_097_hi_hops_netcal-LMTcal-norm_StokesI.uvfits -p ./params/eht-imaging_097_params.csv -c 0 --savepdf
    
For additional details, please read the help document associated in the imaging script: "python eht-imaging_pipeline.py --help".

**Additional References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [Chael, A., Johnson, M., Bouman, K. et al. 2018, ApJ, 857, 23C](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...23C/abstract)
- [Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.2.4)](https://zenodo.org/record/6519440)
- [eht-imaging: https://github.com/achael/eht-imaging](https://github.com/achael/eht-imaging)
