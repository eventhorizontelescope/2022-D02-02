# StarWarps Sgr A* Stokes I Dynamic Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** May 12, 2022

**Primary Reference:** The Event Horizon Telescope Collaboration, et al. 2022, ApJL (Sgr A* Paper III)

**Data Product Code:** [2022-D01-02](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

The pipeline reconstructs a movie from uvfits files simultaneously released in the EHT website (data release ID: 2022-D01-01) using StarWarps (Bouman et al. 2018), implemented in the python-interfaced imaging package eht-imaging (Chael et al. 2018).

To run the pipeline, specify the input uvfits data file, the prior distribution mean image, the file containing the imaging parameters, and the combination to be used. Additional flags control the output, which is only the reconstructed movie as an HDF5 file by default.

Example call:

    python starwarps_pipeline.py -i ../../data/uvfits_besttime/ER6_SGRA_2017_097_lo_hops_netcal-LMTcal_StokesI.uvfits -m ./prior_mean/ring_blur25uas.fits -p ./params/starwarps_params.csv -c 0 --savemp4
    
For additional details, please read the help document associated in the imaging script: "python starwarps_pipeline.py --help".

**Additional References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [Chael, A., Johnson, M., Bouman, K. et al. 2018, ApJ, 857, 23C](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...23C/abstract)
- [Bouman, K., Johnson, M., Dalca, A., et al. 2018, IEEE Transactions on Computational Imaging](https://ui.adsabs.harvard.edu/abs/2017arXiv171101357B/abstract)
- [Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.2.4)](https://zenodo.org/record/6519440)
- [eht-imaging: https://github.com/achael/eht-imaging](https://github.com/achael/eht-imaging)
