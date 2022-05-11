"""
eht-imaging Sgr A* Stokes I Static Imaging Pipeline for EHT observations in April 2017
Authors: The Event Horizon Telescope Collaboration et al.
Date: May 12, 2022
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2022, ApJL (Sgr A* Paper III)
Data Product Code: 2022-D01-02
Brief Description:
The pipeline reconstructs an image from uvfits files simultaneously
released in the EHT website (data release ID: 2022-D01-01) using
the python-interfaced imaging package eht-imaging (Chael et al. 2018).
To run the pipeline, specify the input uvfits data file, the file containing the imaging parameters,
and the combination to be used. Multiple bands can be specified separately.
Additional flags control the output, which is only the reconstructed image as a FITS file by default.
Example call:
python eht-imaging_pipeline.py -i ../../data/uvfits_norm/ER6_SGRA_2017_097_lo_hops_netcal-LMTcal-norm_StokesI.uvfits -i2 ../../data/uvfits_norm/ER6_SGRA_2017_097_hi_hops_netcal-LMTcal-norm_StokesI.uvfits -p ./params/eht-imaging_097_params.csv -c 0 --savepdf
For additional details, please read the help document associated in the
imaging script: "python eht-imaging_pipeline.py --help".
Additional References:
 - EHT Collaboration Data Portal Website:
   https://eventhorizontelescope.org/for-astronomers/data
 - Chael, A., Johnson, M., Bouman, K., et al. 2018, ApJ, 857, 23C
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
import numpy as np
import pandas as pd
import preimcal

#-------------------------------------------------------------------------------
# Load command-line arguments
#-------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Fiducial eht-imaging script for Sgr A*")
parser.add_argument('-i', '--infile',   default="obs.uvfits", help="input UVFITS file")
parser.add_argument('-i2', '--infile2', default="",           help="optional 2nd input file (different band) for imaging")
parser.add_argument('-p', '--params',   default='',           help='input parameters file')
parser.add_argument('-c', '--combo',    default=0,            help='parameter combination number', type=int)
parser.add_argument('-o', '--outfile',  default='out.fits',   help='output FITS image')
parser.add_argument('--savepdf',        default=False,        help='saves image pdf (True or False)?', action='store_true')
args = parser.parse_args()

#-------------------------------------------------------------------------------
# Load parameter combination from the Top Set images obtained from the
# eht-imaging parameter survey
#-------------------------------------------------------------------------------

pset = pd.read_csv(args.params)
pset = pset.iloc[args.combo]

#-------------------------------------------------------------------------------
# Fiducial imaging parameters obtained from the eht-imaging parameter survey
#-------------------------------------------------------------------------------

# Do network calibration
netcal       = False

# Do LMT calibration
LMTcal       = False
LMTcal_fwhm  = 60*eh.RADPERUAS

# Zero baseline flux density
zbl          = 1

# Light curve array
lcarr        = None

# Do JCMT calibration
JCMTcal      = False

# Are input data normalized
normalized   = True

# Are input data deblurred
deblurred    = False

# Time averaging interval (s)
tint         = 60

# Systematic noise budget
syserr       = pset.syserr

# Gains solution interval
solint       = 0

# Refractive scattering noise type
ref_type     = pset.ref_type

# Refractive scattering noise scaling
ref_scale    = pset.ref_scale

# Do data deblurring
deblurr      = pset.deblurr

# PSD noise model (see Paper III, Eq. 2)
psd_noise    = pset.psd_noise
psd_a        = pset.psd_a
psd_u0       = pset.psd_u0
psd_b        = pset.psd_b
psd_c        = pset.psd_c

# Gains tolerance during self-calibration
gaintol      = (0.02, 0.2)

# Do self-calibration
selfcal      = False

# Field of view (uas)
fov          = 150*eh.RADPERUAS

# Number of pixels to reconstruct
npixels      = 80

# Second moment Gaussian parameters
xmaj         = 60*eh.RADPERUAS  # major axis (uas)
xmin         = 60*eh.RADPERUAS  # minor axis (uas)
xpa          = 0*eh.RADPERUAS   # position angle (deg)

# Prior/initialization image size 
fwhm         = pset.fwhm*eh.RADPERUAS  # Prior size (uas)

# Number of imaging iterations
niter_static = 5

# Initialization image blurring fraction
blurfrac     = 1

# Number of convergence iterations
maxit        = 500

# Data products weighting
amp          = pset.amp     # Visibility amplitudes
cphase       = 1            # Closure phases
logcamp      = 1            # Log closure amplitudes

# Regularizers weighting
simple       = pset.simple  # Maximum-Entropy
tv           = pset.tv      # Total variation
tv2          = pset.tv2     # Total squared variation
flux         = 1            # Compact flux constraint

# Epsilon parameter for total variation
epsilon_tv   = 1e-10

# Convergence criterion
stop         = 1e-6

# Type of Fourier transform ('direct', 'nfft', or 'fast')
ttype        = 'nfft'

# Regularization terms dictionary
reg_term     = {'simple' : simple,
                'tv'     : tv,
                'tv2'    : tv2,
                'flux'   : flux}

# Data products dictionary
data_term    = {'amp'    : amp,
                'cphase' : cphase,
                'logcamp': logcamp}

#-------------------------------------------------------------------------------
# Prepare the data
#-------------------------------------------------------------------------------

obsfile = args.infile  # input data
nproc   = -1           # parallel processes (-1 = all available)
    
# Load the input uvfits file
obs1 = eh.obsdata.load_uvfits(obsfile)

# Prepare the data
inputset = [obs1, normalized, deblurred, lcarr, nproc, LMTcal, 
            LMTcal_fwhm, JCMTcal, tint, syserr, ref_type, ref_scale,
            deblurr, psd_noise, psd_a, psd_u0, psd_b, psd_c]

print('\n Preparing Input Data')
obs1 = preimcal.preim_pipeline(*inputset)
obs1 = obs1.switch_polrep('stokes')

# If two uvfits files are passed as input (e.g., high and low band) then use
# both datasets, but do not form closure quantities between the two datasets
if args.infile2 != '':
    obsfile2 = args.infile2

    # Load the second input uvfits file
    obs2 = eh.obsdata.load_uvfits(obsfile2)

    # Prepare the data
    inputset = [obs2, normalized, deblurred, lcarr, nproc, LMTcal, 
                LMTcal_fwhm, JCMTcal, tint, syserr, ref_type, ref_scale,
                deblurr, psd_noise, psd_a, psd_u0, psd_b, psd_c]

    print('\n Preparing Input Data 2')
    obs2 = preimcal.preim_pipeline(*inputset)
    obs2 = obs2.switch_polrep('stokes')

    # add a slight offset to avoid mixed closure products
    obs2.data['time'] += 0.00002718

    # concatenate the observations into a single observation object
    print('\n Merging Input Data 1 and 2')
    obs = obs1.copy()
    obs.data = np.concatenate([obs1.data, obs2.data])
    
else:
    obs = obs1.copy()

#-------------------------------------------------------------------------------
# Set up an initial and prior image
#-------------------------------------------------------------------------------

res = obs.res()  # The nominal array resolution: 1/(longest baseline)

# Make a disk prior image for maximum entropy regularization
# This is also the initial image
diskprior = eh.image.make_square(obs, npixels, fov)
diskprior = diskprior.add_tophat(zbl, fwhm/2.)
diskprior = diskprior.blur_circ(res)

#-------------------------------------------------------------------------------
# Define helper function
#-------------------------------------------------------------------------------

# Repeat imaging with blurring to assure good convergence
def converge(major=niter_static, blur_frac=blurfrac):
    for repeat in range(major):
        imgr.make_image_I(show_updates=False)
        init = imgr.out_last().blur_circ(blur_frac*res)
        imgr.init_next = init

#-------------------------------------------------------------------------------
# Reconstruct an image
#-------------------------------------------------------------------------------

print("\n Imaging with visibility amplitudes and closure quantities...")

# Set up imager
imgr = eh.imager.Imager(obs, diskprior, prior_im=diskprior,
                        flux=zbl, 
                        data_term=data_term, reg_term=reg_term,
                        maxit=maxit, norm_reg=True,
                        epsilon_tv=epsilon_tv,
                        ttype=ttype,
                        major=xmaj, 
                        minor=xmin, 
                        PA=xpa,
                        stop=stop)
                        
# Imaging
converge()

#-------------------------------------------------------------------------------
# Output the results
#-------------------------------------------------------------------------------

# Save the reconstructed image
im_out = imgr.out_last().copy()
im_out.save_fits(args.outfile)

# Optionally save a pdf of the final image
if args.savepdf:
    pdfout = os.path.splitext(args.outfile)[0] + '.pdf'
    im_out.display(cbar_unit=['Tb'], label_type='scale', export_pdf=pdfout)
    
