#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = '2.0'

import os
import argparse
# combining bands
from uv_comb import uvf_combine


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FUNCTIONS
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
def merge_bands(inpath, obsfile):
    """
    Merge LO and HI band data sets
    """

    if '_LO_' in obsfile:
        obsfile_lo = obsfile
        obsfile_hi = obsfile.replace('_LO_', '_HI_')
    if '_HI_' in obsfile:
        obsfile_lo = obsfile.replace('_HI_', '_LO_')
        obsfile_hi = obsfile

    if '_lo_' in obsfile:
        obsfile_lo = obsfile
        obsfile_hi = obsfile.replace('_lo_', '_hi_')
    if '_hi_' in obsfile:
        obsfile_lo = obsfile.replace('_hi_', '_lo_')
        obsfile_hi = obsfile

    if '_LO_' in obsfile_lo:
        outfile = obsfile_lo.replace('_LO_', '_lo+hi_')
    if '_lo_' in obsfile_lo:
        outfile = obsfile_lo.replace('_lo_', '_lo+hi_')

    uvf_combine([inpath+obsfile_lo+'.uvfits', inpath+obsfile_hi+'.uvfits'], outp=inpath+outfile+'.uvfits')
    return outfile

def run_difmap(obsfile, inpath, outpath, mask_size, ALMA_weight, uvbin):
    """
    Run Difmap imaging
    """
    # fixed parameters
    rms_target = 0.0001
    source_flux = 2.4
    clean_niter = 100
    clean_gain = 0.2
    uvpower = -1

    if os.path.exists(inpath+obsfile+'.uvfits'):
        print('\n MERGED DATA PROPERLY CREATED! RUNNING DIFMAP\n')
        # run scripted Difmap
        cmd = 'echo @EHT_FitClean %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s | difmap'%(inpath, obsfile,
                                                                                        outpath, mask_size,
                                                                                        rms_target, source_flux,
                                                                                        clean_niter, clean_gain,
                                                                                        ALMA_weight, uvbin, uvpower)
        os.system(cmd)
        os.system('rm difmap.log*')
    else:
        print('\n NO MERGED DATA! STOP DIFMAP\n')


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# IMAGING
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i','--inpath',metavar='Absolute path of input data',type=str,default=None, help='Absolute path of input data')
parser.add_argument('-o','--outpath',metavar='Absolute path to save outputs',type=str,default=None, help='Absolute path to save outputs')

parser.add_argument('--uvfits', help='Input UVFITS file name')
parser.add_argument('--mask_size', help='Mask size')
parser.add_argument('--ALMA_weight', help='ALMA weight for selfcal')
parser.add_argument('--uvbin', help='UV weight, number of bins')

inpath = parser.parse_args().inpath
outpath = parser.parse_args().outpath
obsfile = parser.parse_args().uvfits
mask_size = parser.parse_args().mask_size
ALMA_weight = parser.parse_args().ALMA_weight
uvbin = parser.parse_args().uvbin

if inpath[-1] != '/':
    inpath += '/'
if outpath[-1] != '/':
    outpath += '/'

if not os.path.exists(outpath):
    os.mkdir(outpath)

try:
    # merge LO and HI data sets
    obs_merged = merge_bands(inpath, obsfile)
except:
    print('\n MERGED DATA NOT PROPERLY CREATED!! \n')

# Imaging
run_difmap(obs_merged, inpath, outpath, mask_size, ALMA_weight, uvbin)
