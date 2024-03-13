# First Sgr A* EHT Results: Imaging Pipelines

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** May 12, 2022

**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2022c, ApJL, 930, L14 (Sgr A* Paper III)](https://ui.adsabs.harvard.edu/link_gateway/2022ApJ...930L..14E/doi:10.3847/2041-8213/ac6429)

**Data Product Code:** [2022-D02-01](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

We release three imaging pipelines (DIFMAP, eht-imaging, SMILI, and Starwarps)
used in the parameter survey of Sgr A* Paper III Section 7 and later. All
imaging pipelines create images from calibrated uvfits files (see Sgr A*
Paper II) simultaneously released (data release ID: [2022-D02-01](https://eventhorizontelescope.org/for-astronomers/data)).
For more detailed instructions, please see the README file in the
sub-directory for each pipeline.

We note that, as described in
[2022-D02-01](https://eventhorizontelescope.org/for-astronomers/data),
released visibility data sets have only Stokes *I*, which are slightly
different from data sets used in Paper IV that have dual polarization
at Stokes `RR` and `LL`. This slight difference in released data sets
will provide no net changes in DIFMAP and eht-imaging pipelines, while
it will change self-calibration procedures slightly for SMILI
calibrating `R` and `L` gains separately for the latter dual
polarization data sets. We confirm that reconstructed images are
consistent with images presented in Paper IV on fiducial parameters
(see Paper IV Section 6) for all three pipelines, and will not affect
our conclusions in the M87 publications (Paper I, II, III, IV, V and
VI).

**Notes:**

These data files only include Stokes *I* visibilities, while the
published results used data files from the full SR1 release which
include Stokes `RR` and `LL`. The slight difference in the underlying
data from the conversion to Stokes *I* in the single-precision
`*.uvfits` files, as well as differences in the python dependencies
used by `eht-imaging` and `SMILI` (e.g. `numpy`, `astropy`, `scipy`),
may slightly affect the final image.

**References:**
- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [The Event Horizon Telescope Collaboration, et al. 2022a, ApJL, 930, L12 (Sgr A* Paper I)](https://doi.org/10.3847/2041-8213/ac6674)
- [The Event Horizon Telescope Collaboration, et al. 2022b, ApJL, 930, L13 (Sgr A* Paper II)](https://doi.org/10.3847/2041-8213/ac6675)
- [The Event Horizon Telescope Collaboration, et al. 2022c, ApJL, 930, L14 (Sgr A* Paper III)](https://doi.org/10.3847/2041-8213/ac6429)
- [The Event Horizon Telescope Collaboration, et al. 2022d, ApJL, 930, L15 (Sgr A* Paper IV)](https://doi.org/10.3847/2041-8213/ac6736)
- [The Event Horizon Telescope Collaboration, et al. 2022e, ApJL, 930, L16 (Sgr A* Paper V)](https://doi.org/10.3847/2041-8213/ac6672)
- [The Event Horizon Telescope Collaboration, et al. 2022f, ApJL, 930, L17 (Sgr A* Paper VI)](https://doi.org/10.3847/2041-8213/ac6756)
