"""
StarWarps Sgr A* Stokes I Dynamic Imaging Pipeline for EHT observations in April 2017
Authors: The Event Horizon Telescope Collaboration et al.
Date: May 12, 2022
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2022, ApJL (Sgr A* Paper III)
Data Product Code: 2022-D01-02
Brief Description:
The pipeline reconstructs a movie from uvfits files simultaneously
released in the EHT website (data release ID: 2022-D01-01) using StarWarps (Bouman et al. 2018),
implemented in the python-interfaced imaging package eht-imaging (Chael et al. 2018).
To run the pipeline, specify the input uvfits data file, the prior distribution mean image,
the file containing the imaging parameters, and the combination to be used.
Additional flags control the output, which is only the reconstructed movie as an HDF5 file by default.
Example call:
python starwarps_pipeline.py -i ../../data/uvfits_besttime/ER6_SGRA_2017_097_lo_hops_netcal-LMTcal_StokesI.uvfits -m ./prior_mean/ring_blur25uas.fits -p ./params/starwarps_params.csv -c 0 --savemp4
For additional details, please read the help document associated in the
imaging script: "python starwarps_pipeline.py --help".
Additional References:
 - EHT Collaboration Data Portal Website:
   https://eventhorizontelescope.org/for-astronomers/data
 - Chael, A., Johnson, M., Bouman, K., et al. 2018, ApJ, 857, 23C
 - Bouman, K., Johnson, M., Dalca, A., et al. 2018, IEEE Transactions on Computational Imaging
 - Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.2.4)
 - eht-imaging: https://github.com/achael/eht-imaging
"""

#-------------------------------------------------------------------------------
# Author Information
#-------------------------------------------------------------------------------

__author__ = "The Event Horizon Telescope Collaboration et al."
__copyright__ = "Copyright 2022, the Event Horizon Telescope Collaboration et al."
__license__ = "GPL version 3"
__version__ = "1.0"
__date__  = "May 12 2022"

#-------------------------------------------------------------------------------
# Modules
#-------------------------------------------------------------------------------

import matplotlib
matplotlib.use('Agg')

import os
import argparse
import ehtim as eh
from ehtim.imaging import starwarps as sw
import numpy as np
import pandas as pd
import preimcal
import contextlib
import scipy.interpolate as interp

#-------------------------------------------------------------------------------
# Load command-line arguments
#-------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Fiducial StarWarps script for Sgr A*")
parser.add_argument('-i', '--infile',   default="obs.uvfits", help="input UVFITS file")
parser.add_argument('-m', '--mean',     default="mean.fits",  help="input prior mean FITS file")
parser.add_argument('-p', '--params',   default='',           help='input parameters file')
parser.add_argument('-c', '--combo',    default=0,            help='parameter combination number', type=int)
parser.add_argument('-o', '--outfile',  default='out.hdf5',   help='output HDF5 movie')
parser.add_argument('--savemp4',        default=False,        help='saves movie mp4 (True or False)?', action='store_true')
parser.add_argument('--movsum',         default=False,        help='generate movie summary pdf', action='store_true')
args = parser.parse_args()

#-------------------------------------------------------------------------------
# Load parameter combination from the StarWarps parameter survey
#-------------------------------------------------------------------------------

pset = pd.read_csv(args.params)
pset = pset.iloc[args.combo]

#-------------------------------------------------------------------------------
# Dynamic imaging parameters used in the StarWarps parameter survey
#-------------------------------------------------------------------------------

# Do network calibration
netcal            = False

# Do LMT calibration
LMTcal            = False
LMTcal_fwhm       = 60*eh.RADPERUAS

# Zero baseline flux density
zbl               = 2.3

# Light curve array
lcarr             = None

# Do JCMT calibration
JCMTcal           = False

# Are input data normalized
normalized        = False

# Are input data deblurred
deblurred         = False

# Time averaging interval
tint              = 60.

# Systematic noise budget
syserr            = 0.02

# Gains solution interval
solint            = 0

# Refractive scattering noise type
ref_type          = pset.ref_type

# Refractive scattering noise scaling
ref_scale         = pset.ref_scale

# Do data deblurring
deblurr           = pset.deblurr

# PSD noise model (see Paper III, Eq. 2)
psd_noise         = False
psd_a             = 0.02
psd_u0            = 2
psd_b             = 2
psd_c             = 2

# Gains tolerance during self-calibration
gaintol           = (0.02, 0.2)

# Do self-calibration
selfcal           = False

# Field of view (uas)
fov               = 150*eh.RADPERUAS

# Number of pixels to reconstruct
npixels           = 40

# Re-scale prior total flux
rescale           = True

# Flow field parameterization
warp_method       = 'phase'

# Number of re-linearization iterations
numLinIters       = 5

# Set prior for every frame (True) or just the first frame (False)
interiorPriors    = True

# Reassign images at each linearization step
reassign_apxImgs  = False

# Temporal regularization term (larger values allow more variability)
variance_img_diff = pset.variance_img_diff

# Covariance smoothness term (larger values are smoother)
powerDropoff      = pset.powerDropoff

# Constraint on fraction of pixels > 0
covfrac           = pset.covfrac

# Data products weighting
amp               = pset.amp      # Visibility amplitudes
logcamp           = pset.logcamp  # Log closure amplitudes
bs                = 10            # Bispectrum
flux              = 1e4           # Compact flux constraint

# Type of Fourier transform ('direct', 'nfft', or 'fast')
ttype             = 'nfft'

# Data products dictionary
data_term         = {'amp'    : amp,
                     'logcamp': logcamp,
                     'bs'     : bs,
                     'flux'   : flux}

#-------------------------------------------------------------------------------
# Prepare the data
#-------------------------------------------------------------------------------

obsfile = args.infile  # input data
nproc   = -1           # parallel processes (-1 = all available)
    
# Load the input uvfits file
obs = eh.obsdata.load_uvfits(obsfile)

# Prepare the data
inputset = [obs, normalized, deblurred, lcarr, nproc, LMTcal, 
            LMTcal_fwhm, JCMTcal, tint, syserr, ref_type, ref_scale,
            deblurr, psd_noise, psd_a, psd_u0, psd_b, psd_c]

print('\n Preparing Input Data')
obs = preimcal.preim_pipeline(*inputset)
obs = obs.switch_polrep('stokes')

#-------------------------------------------------------------------------------
# Set up a prior mean image
#-------------------------------------------------------------------------------

res   = obs.res()  # The nominal array resolution: 1/(longest baseline)
prior = eh.image.load_fits(args.mean)  # Input prior distribution mean image

# Ensure rest frequency, right-ascension, and declination info consistency
prior.rf  = obs.rf
prior.ra  = obs.ra
prior.dec = obs.dec

# Regrid to correct resolution and pixel number
prior = prior.regrid_image(fov, npixels)

# Rescale total flux
if rescale:
    prior.imvec *= zbl / prior.total_flux()

#-------------------------------------------------------------------------------
# Obtain light curve from the data
#-------------------------------------------------------------------------------

def get_lightcurve(obs):
    obs.add_scans()
    obs_split = obs.split_obs(scan_gather=False)

    # Univariate spline interpolation
    alltimes = np.array([obs_split[j].data['time'][0] for j in range(len(obs_split))])
    with contextlib.redirect_stdout(None):
        allfluxes = np.array(
            [np.median(obs_split[j].flag_uvdist(uv_max=0.01e9).unpack('amp')['amp']) 
             for j in range(len(obs_split))])

    idxsort = np.argsort(alltimes)
    alltimes = alltimes[idxsort]
    allfluxes = allfluxes[idxsort]

    mask = np.isnan(allfluxes)
    maskedtimes = alltimes[~mask]
    maskedfluxes = allfluxes[~mask]

    spl = interp.UnivariateSpline(maskedtimes, maskedfluxes, ext=3)
    spl.set_smoothing_factor(1e-10)
    spl_times = alltimes
    spl_fluxes = spl(spl_times)
    flux_list = list(spl_fluxes)

    return flux_list

#-------------------------------------------------------------------------------
# Reconstruct a movie
#-------------------------------------------------------------------------------

# Initalize lists
imCov, meanImg, initImg = [], [], []

meanImg.append(prior) 
initImg.append(prior)

# Set covariance
imCov.append(sw.gaussImgCovariance_2(meanImg[0],
                                     powerDropoff=powerDropoff,
                                     frac=covfrac))
                                     
# Make the covariance matrix that says how much variation there should be
# between frames in time
noiseCov_img = np.eye(npixels**2) * variance_img_diff


# Initialize the flowbasis and get the initTheta which says how to
# specify no motion for the specified flow basis
init_x, init_y, flowbasis_x, flowbasis_y, initTheta = \
    sw.affineMotionBasis_noTranslation(meanImg[0])

# Split observation into specified time intervals
obs_list = obs.split_obs(t_gather=tint)
                   
# Get data light curve
flux_list = get_lightcurve(obs)

# Iterate
frames, expVal_t_t, expVal_tm1_t, loglikelihood, apxImgs = \
    sw.computeSuffStatistics(meanImg, imCov,
                obs_list, noiseCov_img,
                initTheta, init_x, init_y,
                flowbasis_x, flowbasis_y, initTheta,
                method=warp_method,
                measurement=data_term,
                init_images=initImg, 
                lightcurve=flux_list,
                interiorPriors=interiorPriors,
                numLinIters=numLinIters,
                compute_expVal_tm1_t=False)

# Asign the time corresponding to each frame
for i in range(len(frames)):
    frames[i].time = obs_list[i].data['time'][0]
    frames[i].mjd  = obs_list[i].mjd

# Merge frames into a movie
movie = eh.movie.merge_im_list(frames)
movie.reset_interp(interp='linear', bounds_error=False)

# Self-calibrate data
if selfcal:
    obs_movie = eh.selfcal(obs, movie, method='both', ttype=ttype,
                           solution_interval=solint, gain_tol=gaintol,
                           processes=nproc)
else:
    obs_movie = obs.copy()

#-------------------------------------------------------------------------------
# Output the results
#-------------------------------------------------------------------------------

# Save the reconstructed movie
movie.save_hdf5(args.outfile)

# Optionally save a mp4 of the final movie
if args.savemp4:
    mp4out = os.path.splitext(args.outfile)[0] + '.mp4'
    movie.export_mp4(out=mp4out, label_time=True, fps=10)

# Optionally generate a summary of the final movie and associated data consistency metrics
if args.movsum:
    # Save an movie summary sheet
    matplotlib.pyplot.close('all')
    outmovsum = os.path.splitext(args.outfile)[0] + '_movsum.pdf'
    eh.imgsum(movie, obs_movie, obs, outmovsum, processes=nproc)
    
