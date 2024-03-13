#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""SMILI Sgr A* Stokes I Imaging Pipeline for EHT observations in April 2017

Authors: The Event Horizon Telescope Collaboration et al.
Date: May 10, 2022
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2022c, ApJL, 930, L14 (Sgr A* Paper III)
Data Product Code: 2022-D02-01

Brief Description:
The pipeline reconstructs an image from a specified uvfits file simultaneously
released in the EHT website (data release ID: 2022-D02-01) using a
python-interfaced imaging package SMILI (Akiyama et al. 2017a,b) version 0.2.0
(Moriyama et al. 2022).This script assumes to load calibrated data sets at Stokes
I (see M87 Paper III). As described in M87 Paper IV Section 6.2.3, the SMILI
pipeline also uses the eht-imaging library (Chael et al. 2016, 2018)
version 1.2.3 (Chael et al. 2019) solely for pre-imaging calibrations including
time averaging, LMT calibration to use the input visibility data sets consistent
with the eht-imaging imaging pipeline. The script has been tested in Python 3
installed with Anaconda on Ubuntu 18.04 LTS and MacOS 10.13 & 10.14 (with
macport or homebrew).

Notes:
The pipeline requires the pre-imaging calibration script (preimcal.py) and scattering pickle file
(e.g., obs_scatt_std_gauss_scaled_3599.pickle) at the same directory.

The pipeline will output three files.
 - image (specified with -o option; assumed to be xxxx.fits here)
 - pre-calibrated uvfits file (xxxx.precal.uvfits)
 - self-calibrated uvfits file (xxxx.selfcal.uvfits)
For usage and detail parameters of the pipeline, please
read the help document associated in the imaging script
"python smili_imaging_pipeline.py --help".

References:
 - EHT Collaboration Data Portal Website:
   https://eventhorizontelescope.org/for-astronomers/data
 - The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)
 - The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)
 - The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)
 - The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)
 - The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)
 - The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)
 - Akiyama, K., Ikeda, S., Pleau, M., et al. 2017a, ApJ, 838, 1
 - Akiyama, K., Kuramochi, K., Ikeda, S., et al. 2017b, AJ, 153, 159
 - Moriyama, K., Akiyama, K., Cho, I., et al. 2022, Zenodo (SMILI version 0.2.0)
 - Chael, A., Johnson, M., Narayan, R., et al. 2016, ApJ, 829, 11C
 - Chael, A., Johnson, M., Bouman, K., et al. 2018, ApJ, 857, 23C
 - Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.2.3)
 - Anaconda: https://www.anaconda.com/
 - eht-imaging: https://github.com/achael/eht-imaging
 - SMILI: https://github.com/astrosmili/smili
"""

#-------------------------------------------------------------------------------
# Information of Authors
#-------------------------------------------------------------------------------
__author__ = "The Event Horizon Telescope Collaboration et al."
__copyright__ = "Copyright 2020, the Event Horizon Telescope Collaboration et al. in prep"
__license__ = "GPL version 3"
__version__ = "1.0"
__date__  = "September 16 2020"


#-------------------------------------------------------------------------------
# Modules
#-------------------------------------------------------------------------------
from smili import uvdata,imdata,imaging
import ehtim as eh
import pandas as pd
import numpy as np
import copy
import os
import argparse
import itertools
import imaging_function as imf
import preimcal
import random

#-------------------------------------------------------------------------------
# Loading command line arguments
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i','--input',metavar='Input uvfits file',nargs="+",type=str,required=True,
                    help='Input uvfits file name')
parser.add_argument('-o','--output',metavar='Output fits file',type=str,default=None,
                    help='Output fits file name')
# Scattering
parser.add_argument('--reftype', metavar='Refractive noise type', type=str, default="None",
                    help='Refractive noise table (quarter1, quarter2, dime). If not specified, do not add this budget.')
parser.add_argument("--refscale", metavar="Scale of refractive scattering", type=float, default=1.)
parser.add_argument('--deblurr', action='store_true', default=False,
                    help='If specified, the output will be deblurred.')
# Normalization
parser.add_argument('--is_normalized', action='store_true', default=False,
                    help='If the input uvfits is flux-normalized, please DO NOT FORGET TO SPECIFY THIS OPTION.')
parser.add_argument("--shifttype",metavar='shifttype',type=str,default=None,
                    help='Image shift type (e.g., "peak", "com1", "hough", "nxcorr")')

parser.add_argument('--lmtfwhm',type=float,default=-1,
                    help='FWHM size for the pre-imaging LMT selfcal. If a non-positive value is specified, do not do any preLMT cal.')
# Prior
parser.add_argument("--prior_image", metavar="Prior image structure", type=str, default=None,
                    help="Prior image structure (tophat, gaussian, or input fits file name)")
parser.add_argument("--prior_size",metavar='Prior size',type=float,default=100.,
                    help='Prior size (uas)')
# Time information
parser.add_argument("--tint",metavar='Intgeration Time',type=float,default=-1,
                    help='Integration Time in seconds. Negative value means no averaging, 0 means scan-averaging.')
parser.add_argument("--wint",metavar='Weigthcal Solint',type=float,default=-1,
                    help='Solution Interval of weigthcal in seconds. Negative value means no weightcal.')
# systematic errors and cutoff
parser.add_argument("--syserr",metavar='Fractional systematic error',type=float,default=0.01,
                    help='Fractional systematic error to be added to visibilities in quadrature (percent)')
parser.add_argument('--snrcutoff',type=float,default=-1,
                    help='SNR cutoff for visibilities')
# Imaging regularizers
parser.add_argument("--lambl1",metavar='lambda l1',type=float,default=1.,
                    help='lambda l1')
parser.add_argument("--l1_noise",metavar='Floor of l1 prior',type=float,default=0.01,
                    help='Floor of l1 prior')
parser.add_argument("--lambtv",metavar='lambda tv',type=float,default=1e2,
                    help='lambda tv')
parser.add_argument("--lambtsv",metavar='lambda tsv',type=float,default=1e2,
                    help='lambda tsv')
parser.add_argument("--lambgs",metavar='lambda gs',type=float,default=1e-1,
                    help='lambda gs')
# Initial movie/image
parser.add_argument('--init_im', metavar='Initial image structure (tophat, gaussian, or input fits file name)',type=str, default="tophat",
                    help='Initial image structure (tophat, gaussian, or input fits file name).')
parser.add_argument("--init_im_size",metavar='Initial image size (uas)',type=float,default=20.,
                    help='Initial image size (uas).')

# Dynamical imaging information
#   Time interval between each frame of a movie
parser.add_argument('--time_frame',metavar='Time frame information',type=str,default="static",
                    help='Information for a time frame of a movie (time interval (minutes), or time table name)')
#   Light curve
parser.add_argument('--lc_type',metavar='Light curve information',type=str,default='zbl',
                    help='Light curve information (zbl, uniform total flux value, or table name). Default=zbl (zero baseline).')
parser.add_argument('--extflux', metavar='Extended flux', type=float, default=0.,
                    help='Extended flux (Jy). Light curve (lc) of a movie = original lc - external flux.')
#   Regularizers
parser.add_argument("--wamp",metavar="Weight of visibility amplitudes",type=float,default=1.,
                    help="Weight of visibility amplitudes")
parser.add_argument("--lambrt",metavar="lambda rt",type=float,default=-1,
                    help="lambda rt")
parser.add_argument("--lambri",metavar="lambda ri",type=float,default=-1,
                    help="lambda ri")

# Parameters of the noise modeling
parser.add_argument("--a",metavar="PSD a",type=float,default=-1,
                    help="PSD Noise Param: Scaling factor in Jy. Typically (0, 0.1) Jy.")
parser.add_argument("--u0",metavar="PSD u0",type=float,default=-1,
                    help="PSD Noise Param: Break location of the broken power law in Glambda. Typically (0, 2).")
parser.add_argument("--b",metavar="PSD b",type=float,default=-1,
                    help="PSD Noise Param: power-law index at longer baselines than u0. Typically (0, 6).")
parser.add_argument("--c",metavar="PSD c",type=float,default=-1,
                    help="PSD Noise Param: power-law index at shorter baselines than u0. Typically (1.5, 2.5).")
# Other arguments
parser.add_argument("--pols",metavar='pols',type=str,default="I",
                    help='Polarization (e.g. "I", "RR+LL")')
parser.add_argument("--nproc",metavar="Number of parallel processors. Default=1.",type=int,default=1)
parser.add_argument("--nx",metavar="Number of parallel processors. Default=100",type=int,default=64)
args = parser.parse_args()

# Input uvfits file
inputuvfile_list = args.input
# Output directory
outputfits = args.output
if outputfits is None:
    outputfits = os.path.splitext(inputuvfile)[0]+".fits"
outputhead = os.path.splitext(outputfits)[0]

# Prior size
prior_image = args.prior_image
prior_size = args.prior_size

# Systematic error
syserr = args.syserr

# Lambda L1, TV and TSV
lambl1  = args.lambl1
lambtv  = args.lambtv
lambtsv = args.lambtsv
lambgs = args.lambgs

# PSD noise parameters
a  = args.a
u0  = args.u0
b = args.b
c = args.c

# Dynamical rt and ri regularizer
time_frame = args.time_frame
lambrt      = args.lambrt
lambri      = args.lambri
if time_frame =="static":
    lambrt      = -1
    lambri      = -1

# Scattering
ref_optype = args.reftype
do_ref = ref_optype!="False"
do_deblurr = args.deblurr
refscale = args.refscale

# Time average
do_timeavg = args.tint >= 0
do_weightcal = args.wint >= 0
is_normalized = args.is_normalized
tint = args.tint
wint = args.wint
# LMT calibration
lmtfwhm = args.lmtfwhm
do_LMTcal = lmtfwhm > 0

# SNR curoff
snrcutoff = args.snrcutoff
# Reorder station number for closure quantites using median SNR (False or "snr")
reorder_st = True

# Number of processors
nproc = args.nproc
os.environ['OMP_NUM_THREADS'] = '%d'%(nproc)

# polarization
pols = args.pols.split("+")

# shift type
shifttype = args.shifttype

# light curve or totalflux information
lc_type = args.lc_type
extflux = args.extflux

#-------------------------------------------------------------------------------
# Other tuning parameters fixed in the paper.
#-------------------------------------------------------------------------------
# Number of iterative imaging inside a single selfcal cycle
Nweig = 8

# Number of closure imaging attempted in this script (must be < Nself)
Nflcl = 2

# Number of selfcals
Nself = 2

# Imaging pixel size
dx_uas = 2

# Number of pixels
nx = args.nx

#-------------------------------------------------------------------------------
# Set default imaging parameters
#-------------------------------------------------------------------------------

# Set a prior for L1
l1_prior = imf.set_l1prior(prior_image,prior_size,dx_uas,nx)
# Default imaging parameters
imprm_init={}

# Iteration numbers
imprm_init["niter"] = 5000
# Regularization parameters
#   Flat prior TSV
if lambtsv > 0:
    imprm_init["tsv_lambda"] = lambtsv
else:
    imprm_init["tsv_lambda"] = -1

#   Flat prior TV
if lambtv > 0:
    imprm_init["tv_lambda"] = lambtv
else:
    imprm_init["tv_lambda"] = -1

#   Weighted L1
if lambl1 > 0:
    imprm_init["l1_lambda"] = lambl1
    imprm_init["l1_prior"] = l1_prior
    imprm_init["l1_noise"] = args.l1_noise

else:
    imprm_init["l1_lambda"] = -1
    imprm_init["gs_prior"] = l1_prior

#   Maximum Entropy Methods: not used in SMILI pipeline
#     GS-MEM
imprm_init["gs_lambda"] = lambgs
#     KL-MEM
imprm_init["kl_lambda"] = -1

# Dynamical regularizer
imprm_init["rt_lambda"] = lambrt
imprm_init["ri_lambda"] = lambri
movietime = 2. # total time of output movie [sec]

imprm_init["nprint"] = 500
imprm_init["inorm"] = 1.

# Initial image (movie)
init_im      = args.init_im
init_im_size = args.init_im_size

#-------------------------------------------------------------------------------
# Step 1-1: Pre-calibration
#   To make sure that we use visibility data sets pre-calibrated consistently with
#   two other pipelines, this pipeline also uses eht-imaging library.
#-------------------------------------------------------------------------------

# Load uvfits data (optional precalibration (LMTcal, JCMTcal, and reweight time))
obs_list_precal1=[]
uvfits_list_precal1 = []
for inputuvfile in inputuvfile_list:
    obs = eh.obsdata.load_uvfits(inputuvfile)
    obs = preimcal.preim_pipeline(obs, is_normalized=is_normalized,
                    do_LMTcal=do_LMTcal, LMTcal_fwhm=lmtfwhm, tint=tint,
                    ref_optype=False, do_deblurr=False,do_psd_noise=False)

    obs_list_precal1.append(obs)

    uvfile_tmp="./tmp_precal1%08d.uvfits"%(random.randint(0,1e+8))
    obs.save_uvfits(uvfile_tmp)
    uvfits = uvdata.UVFITS(uvfile_tmp)
    uvfits_list_precal1.append(uvfits)
    # set movie frame
    if len(obs_list_precal1)==1:
        if time_frame=="static":
            initmovie = imf.set_initialmovie_pipeline(init_im,init_im_size,time_frame,inputuvfile_list[0],dx_uas,nx)
        else:
            initmovie = imf.set_initialmovie_pipeline(init_im,init_im_size,time_frame,uvfile_tmp,dx_uas,nx, "non-uniform")
            imprm_init["rt_prior"] = l1_prior
        Nt = initmovie.Nt
    os.system("rm %s"%(uvfile_tmp))
# Set light curve after time average
lctable = imf.set_lightcurve_movie(initmovie,lc_type,extflux,inputuvfile_list)

imprm_init["lc_lambda"] = 1
imprm_init["lc_tgterror"] = 0.01
imprm_init["lc_array"] = lctable["flux"]
# Rescaling of coefficients of dynamical regularizer with input light curve
imprm_init["lc_normalize"] = "dynamic"

#---------------------------------------------------------------------------
# Step 2: Generating Data Tables for Imaging
#---------------------------------------------------------------------------
# Create the full complex visibility table
uvfits_list_precal2 = []

for obs in obs_list_precal1:
    obs = preimcal.preim_pipeline(obs, is_normalized=is_normalized, syserr=syserr,
                            ref_optype=ref_optype, ref_scale=refscale, do_deblurr=do_deblurr,
                            a=a, u0=u0, b=b, c=c)
    uvfile_tmp="./tmp_precal1%08d.uvfits"%(random.randint(0,1e+8))
    obs.save_uvfits(uvfile_tmp)
    uvfits = uvdata.UVFITS(uvfile_tmp)
    os.system("rm %s"%(uvfile_tmp))
    uvfits_list_precal2.append(uvfits)

# Construct visibilty tables
vtable_list, btable_list, ctable_list = [], [], []
for pol,uvfits in itertools.product(pols,uvfits_list_precal2):
    vtable = uvfits.select_stokes(pol).make_vistable()
    if snrcutoff > 0:
        vtable = vtable.snrcutoff(threshold=snrcutoff).sort_values(by="utc").reset_index(drop=True)
    vtable_list.append(vtable)
    btable_list.append(vtable.make_bstable_min(redundant=[["AA","AP"],["JC","SM"]],reorder_st=reorder_st))
    ctable_list.append(vtable.make_catable(redundant=[["AA","AP"],["JC","SM"]],reorder_st=reorder_st))
vtable = initmovie.set_frmidx(pd.concat(vtable_list)).sort_values(by="utc").reset_index(drop=True)
btable = initmovie.set_frmidx(pd.concat(btable_list)).sort_values(by="utc").reset_index(drop=True)
ctable = initmovie.set_frmidx(pd.concat(ctable_list)).sort_values(by="utc").reset_index(drop=True)
del vtable_list, btable_list, ctable_list

# Check the number of data sets
if len(vtable.amp) == 0:
    raise ValueError("No data points exist in the input visibility data set.")

# Flag data outside of the input time table
idx = np.where(vtable.frmidx.values>=0)
vtable = copy.deepcopy(vtable.loc[idx[0], :])
idx = np.where(btable.frmidx.values>=0)
btable = copy.deepcopy(btable.loc[idx[0], :])
idx = np.where(ctable.frmidx.values>=0)
ctable = copy.deepcopy(ctable.loc[idx[0], :])
del idx

# Check the number of data sets
if len(btable.amp) == 0:
    btable = None
if len(ctable.amp) == 0:
    ctable = None

# fit beam
beamprm = vtable.fit_beam(imdata.IMFITS(dx=dx_uas, nx=nx, angunit="uas"))
cbeamprm = beamprm.copy()
cbeamprm["majsize"] = (beamprm["majsize"]+beamprm["minsize"])/2.
cbeamprm["minsize"] = cbeamprm["majsize"]


#---------------------------------------------------------------------------
# Step 3: Imaging
#---------------------------------------------------------------------------

# Loop for self-calibration
for iself in range(Nself):

    # Iterative Imaging:
    #   At each iteration of self-calibration, imaging iterations are attempted
    #   for Nweig(=5) times.
    for iweig in range(Nweig):
        # Imaging parameters: copied from the default parameter sets
        imprm = copy.deepcopy(imprm_init)

        # Set the initial movie
        for it in range(Nt):
            initmovie.images[it].set_beam(**cbeamprm)

        # Data to be used in default
        #   Closure phases/amplitudes are used always to avoid over-fitting to
        #   full complex visibilities
        imprm["bstable"]=btable
        imprm["catable"]=ctable

        # Use amplitudes or full complex visibilities:
        #   Parameters depending on self-cal stages
        if iself < Nflcl:
            # The first half for self-calibration: using closure quantities

            # Get Amplitude Data sets
            atable = copy.deepcopy(vtable).debias_amp()

            # Flag the intra-site to avoid conflict with with the total-flux-density regularization.
            atable = atable.query("uvdist > 0.1e9").sort_values(by="utc").reset_index(drop=True)

            # Adding a fractional error to the data
            #  For LMT (adding 30% error)
            idx = atable["st1name"] == "LM"
            idx|= atable["st2name"] == "LM"
            atable.loc[idx,"sigma"] = np.sqrt(atable.loc[idx,"sigma"]**2+(atable.loc[idx,"amp"]*0.3)**2)
            #  For other stations (adding an error specified with amp_err)
            idx = idx == False
            atable.loc[idx,"sigma"] = np.sqrt(atable.loc[idx,"sigma"]**2+(atable.loc[idx,"amp"]*0.05)**2)

            # Use amplitudes
            imprm["w_amp"]  = args.wamp
            imprm["amptable"]=atable
        else:
            # The latter half for self-calibration: using full complex
            # visibilities. Note that we still add closure quantities
            # to prevent over-fitting to full complex visibilities
            vtable_in = vtable
            imprm["vistable"] = vtable_in

        # Imaging
        outmovie = initmovie.copy()
        outmovie = imaging.lbfgs3d.imaging3d(outmovie, **imprm)
        initmovie = outmovie.copy()
        initmovie = imf.edit_movie(initmovie,imprm,cbeamprm,shifttype).copy()

    # If iself==Nself, self-calibration won't be performed and exit the script.
    if Nself==1:
        iuv=0
        for uvfits in uvfits_list_precal2:
            if len(uvfits_list_precal2)>1:
                uvfits.to_uvfits(outputhead+"%d.uvfits"%(iuv))
            else:
                uvfits.to_uvfits(outputhead+".uvfits")
            iuv+=1
    if iself == Nself-1:
        break

    # Self-calibration
    uvfits_list_mod=[]
    i=0
    for uvfits in uvfits_list_precal2:
        # 1. Edit images
        #   Before selfcal, we bring the center of the mass to the image center
        #   and also normalize the image with the target total flux density.
        #outimage = outimage.peakshift()
        #   Normalize the total flux density
        #   *** Thanks to the total-flux-density regularization, this will provide only a
        #   *** slight change in the total flux density by a few pecent or even < 1%
        for it in range(Nt):
            outmovie.images[it].data *= lctable["flux"][it]/outmovie.images[it].totalflux()
            outmovie.images[it].set_beam(**cbeamprm)

        # 2. Do selfcal
        #   This is a regularized selfcal function using Gaussian priors
        #   for gain amplitudes and phases. std_amp is the standard deviation of
        #   the Gaussian prior for amplitude gains. std_amp=10000 works as almost
        #   the flat prior, so it is essentially selfcal without any gain constraints.
        output = outmovie.selfcal([uvfits], std_amp=1.0)
        # selfcal uvfits should be the precal-1 data (i.e., LMTcal+JCMTcal+reweight time)
        uvfits_precal1 = copy.deepcopy(uvfits_list_precal1[i])
        # modify the slight difference of time tables
        uvfits_precal1.visdata.coord = uvfits.visdata.coord
        splituvfits = outmovie.split_uvfits([uvfits_precal1])
        uvfits_mod = imdata.movie.apply_cltable(splituvfits[0], splituvfits[1], output[2])[0]
        uvfits_list_mod.append(uvfits_mod)

        i+=1

    # Update precal1 uvfits (before adding expected errors)
    uvfits_list_precal1 = copy.deepcopy(uvfits_list_mod)
    del uvfits_list_mod, output, uvfits

    # 5. Get a new visibility table
    # update data with preimaging calibration
    iuv=0
    uvfits_list_precal2 = []
    for uvfits in uvfits_list_precal1:
        uvfile_tmp="./tmp_precal1%08d.uvfits"%(random.randint(0,1e+8))
        uvfits.to_uvfits(uvfile_tmp)
        obs = eh.obsdata.load_uvfits(uvfile_tmp) # debug
        os.system("rm %s"%(uvfile_tmp))

        obs = preimcal.preim_pipeline(obs, is_normalized=is_normalized, syserr=syserr,
                                     ref_optype=ref_optype, ref_scale=refscale, do_deblurr=do_deblurr,
                                     a=a, u0=u0, b=b, c=c)
        uvfile_tmp="./tmp_precal1%08d.uvfits"%(random.randint(0,1e+8))
        obs.save_uvfits(uvfile_tmp)
        uvfits = uvdata.UVFITS(uvfile_tmp) # debug
        os.system("rm %s"%(uvfile_tmp))

        uvfits_list_precal2.append(uvfits)

        # save final selfcalibrated data with deblur
        print("Save final selfcal data")
        if iself==Nself-2:
            if len(uvfits_list_precal1)>1:
                uvfits.to_uvfits(outputhead+"%d.uvfits"%(iuv))
            else:
                uvfits.to_uvfits(outputhead+".uvfits")
        iuv+=1

    # Construct visibilty tables
    vtable_list, btable_list, ctable_list = [], [], []
    for pol,uvfits in itertools.product(pols,uvfits_list_precal2):
        vtable = uvfits.select_stokes(pol).make_vistable()
        if snrcutoff > 0:
            vtable = vtable.snrcutoff(threshold=snrcutoff).sort_values(by="utc").reset_index(drop=True)
        vtable_list.append(vtable)
        btable_list.append(vtable.make_bstable_min(redundant=[["AA","AP"],["JC","SM"]],reorder_st=reorder_st))
        ctable_list.append(vtable.make_catable(redundant=[["AA","AP"],["JC","SM"]],reorder_st=reorder_st))
    vtable = initmovie.set_frmidx(pd.concat(vtable_list)).sort_values(by="utc").reset_index(drop=True)
    btable = initmovie.set_frmidx(pd.concat(btable_list)).sort_values(by="utc").reset_index(drop=True)
    ctable = initmovie.set_frmidx(pd.concat(ctable_list)).sort_values(by="utc").reset_index(drop=True)
    del vtable_list, btable_list, ctable_list

    # Flag data outside of the input time table
    idx = np.where(vtable.frmidx.values>=0)
    vtable = copy.deepcopy(vtable.loc[idx[0], :])
    idx = np.where(btable.frmidx.values>=0)
    btable = copy.deepcopy(btable.loc[idx[0], :])
    idx = np.where(ctable.frmidx.values>=0)
    ctable = copy.deepcopy(ctable.loc[idx[0], :])
    del idx

# Save movie
if time_frame !="static":
    print("Save movie")
    outmovie.to_hdf5(outputhead+".hdf5",hdf5type="ehtim")

print("Save averaged image")
outmovie.average().to_fits(outputhead+".fits")
