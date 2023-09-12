# DIFMAP Sgr A* Stokes I Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** May 12, 2022

(should be updated)
**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)

(should be updated)
**Data Product Code:** [2019-D02-01](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**
DIFMAP script for imaging EHT data 

Requires Caltech's DIFMAP software for CLEAN imaging reconstruction, version v2.5n or above, that can be obtained from [ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html](ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html)

To reproduce fiducial EHT DIFMAP images of Sgr A*:

1. Place i) Difmap script EHT_FitClean, ii) python wrapper script difmap_wrapper.py, and iii) bash to execute run_imaging.bash in the same directory

2. Edit run_imaging.bash with your i) path to the input UVFITS file, ii) output path name, iii) input uvfits file name (without extension; assuming the file has extension .uvfits).

3. TopSet imaging parameters including the pre-imaging considerations are in tables under ./hdf5, so the imaging pipeline will load the tables by default. Set param_search_id='all' to run imaging to all TopSet parameters, otherwise specify the parameter id (e.g., param_search_id='0,1,2,3'

4. Use the calling sequence:

    bash run_imaging.bash

The header of the EHT_FitClean file contains further information regarding the choice of script parameters. 

(should be updated)
**References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)](https://doi.org/10.3847/2041-8213/ab0ec7)
- [The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)](https://doi.org/10.3847/2041-8213/ab0c96)
- [The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)](https://doi.org/10.3847/2041-8213/ab0c57)
- [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)
- [The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)](https://doi.org/10.3847/2041-8213/ab0f43)
- [The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)](https://doi.org/10.3847/2041-8213/ab1141)
- [Shepherd, M., Pearson, T., & Taylor, G. B., 1994, BAAS, 26, 987S](https://ui.adsabs.harvard.edu/abs/1994BAAS...26..987S/abstract)
- [Shepherd, M., Pearson, T., & Taylor, G. B., 1995, BAAS, 27, 903S](https://ui.adsabs.harvard.edu/abs/1995BAAS...27..903S/abstract)
- [DIFMAP: ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html](ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html)
