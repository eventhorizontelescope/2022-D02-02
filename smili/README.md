# SMILI Sgr A* Stokes I Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** May 12, 2022

**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2022c, ApJL, 930, L14 (Sgr A* Paper III)](https://ui.adsabs.harvard.edu/link_gateway/2022ApJ...930L..14E/doi:10.3847/2041-8213/ac6429)

**Data Product Code:** [2022-D02-01](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

The pipeline reconstructs an image from a specified uvfits file, which was released simultaneously with this pipeline on the EHT website (data release ID: [2022-D02-01](https://eventhorizontelescope.org/for-astronomers/data)) using a python-interfaced imaging package SMILI (Akiyama et al. 2017a,b) version 0.2.0 (Moriyama et al. 2022).

This script is designed to load calibrated data sets at Stokes I (see Sgr A* Paper II). As detailed in Sgr A* Paper III Appendix D.3, the SMILI pipeline utilizes the eht-imaging library (Chael et al. 2016, 2018) version 1.2.3 (Chael et al. 2019) solely for pre-imaging calibrations including time averaging, LMT calibration to use the input visibility data sets consistent with the eht-imaging imaging pipeline. The script has been tested in Python 3 installed with Anaconda on  Ubuntu 18.04 LTS and MacOS 10.13 & 10.14, using either macport or homebrew.

The pipeline will output three files.

- image (specified with -o option; assumed to be xxxx.fits here)
- pre-calibrated uvfits file (xxxx.precal.uvfits)
- self-calibrated uvfits file (xxxx.selfcal.uvfits)

For usage and detail parameters of the pipeline, please read the help document associated in the imaging script "*python smili_pipeline.py --help*".

We also include an example driver *example_driver.py* to run *smili_imaging_pipeline.py* on all released data sets. You can run it by "*python example_driver.py --uvfitsdir < xxxx/uvfits > --nproc < number of the processors >*" For details, please have a look at the help document by "*python example_driver.py --help*".

**Additional References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [Akiyama, K., Ikeda, S., Pleau, M. et al. 2017a, ApJ, 838, 1](https://ui.adsabs.harvard.edu/#abs/2017ApJ...838....1A)
- [Akiyama, K., Kuramochi, K., Ikeda, S. et al. 2017b, AJ, 153, 159](https://ui.adsabs.harvard.edu/#abs/2017AJ....153..159A)
- [Akiyama, K., Tazaki, F., Moriyama, K. et al. 2019, Zenodo (SMILI version 0.0.0)](https://zenodo.org/record/2616725)
- [Chael, A., Johnson, M., Narayan, R. et al. 2016, ApJ, 829, 11C](https://ui.adsabs.harvard.edu/abs/2016ApJ...829...11C/abstract)
- [Chael, A., Johnson, M., Bouman, K. et al. 2018, ApJ, 857, 23C](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...23C/abstract)
- [Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.1.0)](https://zenodo.org/record/2614016)
- [Anaconda: https://www.anaconda.com/](https://www.anaconda.com/)
- [eht-imaging: https://github.com/achael/eht-imaging](https://github.com/achael/eht-imaging)
- [SMILI: https://github.com/astrosmili/smili](https://github.com/astrosmili/smili)
