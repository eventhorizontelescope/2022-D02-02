#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = '2.0'

import os, sys
import tempfile
import subprocess
import pandas as pd
import numpy as np
import time
import argparse
import multiprocessing
import paramsurvey.params

import preimcal
from uv_comb import uvf_combine


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# FUNCTIONS
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

class HiddenPrints:
    """
    Suppresses printing from the loop function
    """

    def __init__(self, hideprint=False):
        self.hideprint=hideprint

    def __enter__(self):
        if self.hideprint:
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
        else:
            pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.hideprint:
            sys.stdout.close()
            sys.stdout = self._original_stdout
        else:
            pass

with HiddenPrints(hideprint=False):
    import ehtim as eh


class FitCleanSurvey:
    """
    Defines & executes a regular CLEAN parameter survey
    """

    def __init__(self, paramset, inpath, outpath):

        # create attributes from input
        for param in paramset:
            setattr(self, param, paramset[param])

        self.inpath = inpath
        self.outpath = outpath
        if not self.normalized: self.outpath = self.outpath.replace('_normalized', '')
        if not os.path.exists(self.outpath): os.makedirs(self.outpath)

        # make tmp dir
        self.cwd = os.getcwd() + '/'
        self.tmpdirobj = tempfile.TemporaryDirectory()
        self.tmpdir = self.tmpdirobj.name + '/'

        # pre-imaging attributes
        self.input_obs = None
        self.lcarr = None
        self.preimproc = 1

        # other attributes
        if hasattr(self, 'init'):
            self.init_type = int(self.init.split('_')[0])
            self.init_size = int(self.init.split('_')[1])

        self.obsfile_lo = f'{self.inpath}{self.uvfits}.uvfits'
        self.obsfile_hi = self.obsfile_lo.replace('_LO', '_HI')
        self.obsfile_lo_pre = f'{self.tmpdir}{self.uvfits}_precal.uvfits'
        self.obsfile_hi_pre = self.obsfile_lo_pre.replace('_LO_', '_HI_')

        if self.deblurr: sct = 'dsct'
        else: sct = 'sct'
        name = self.uvfits.split('_LO')[0]
        self.outfile = f'{name}_lo+hi_{sct}'

    def prepare_data(self, obsfile):
        """
        Prepare the data set
        """

        self.input_obs = eh.obsdata.load_uvfits(obsfile)

        inputset = ['input_obs', 'normalized', 'deblurred', 'lcarr',
                    'preimproc', 'lmtcal', 'lmtcal_fwhm', 'jcmtcal',
                    'tint', 'syserr', 'ref_type', 'ref_scale', 'deblurr',
                    'psd_noise', 'psd_a', 'psd_u0', 'psd_b', 'psd_c']

        inputs = [getattr(self, param) for param in inputset]
        obs_precal = preimcal.preim_pipeline(*inputs)

        if '_LO' in obsfile: obs_precal.save_uvfits(self.obsfile_lo_pre)
        elif '_HI' in obsfile: obs_precal.save_uvfits(self.obsfile_hi_pre)

    def merge_bands(self):
        """
        Merge LO and HI band data sets
        """

        uvf_combine([self.obsfile_lo_pre, self.obsfile_hi_pre],
                    outp=f'{self.tmpdir}{self.outfile}.uvfits')

    def run_difmap(self):
        """
        Execute Difmap (overrides run_difmap() from parent class)
        """

        cmd = f'echo @{self.cwd}EHT_FitClean ' \
              f'{self.tmpdir},{self.outfile},{self.outpath},{self.mask_size},{self.rms_target},' \
              f'{self.source_flux},{self.clean_niter},{self.clean_gain},{self.ALMA_weight},' \
              f'{self.LMT_amp_gains},{self.uvbin},{self.uvpower},{self.id},'\
              f'| difmap'
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL)

    def run(self):
        """
        Run the survey
        """

        # check whether images are already computed
        if os.path.exists(f'{self.outpath}{self.outfile}.fits'):
            print('\n IMAGE ALREADY COMPUTED, SKIPPING COMBINATION! \n')

        else:
            # prepare LO and HI data sets
            self.prepare_data(self.obsfile_lo)
            self.prepare_data(self.obsfile_hi)

            # move to tmp direc
            os.chdir(self.tmpdir)

            try:
                # merge LO and HI data sets
                self.merge_bands()
            except:
                print('\n MERGED DATA NOT PROPERLY CREATED!! \n')

            if os.path.exists(f'{self.tmpdir}{self.outfile}.uvfits'):
                print('\n MERGED DATA PROPERLY CREATED! RUNNING DIFMAP\n')
                # run scripted Difmap
                self.run_difmap()
            else:
                # write failed models to log
                with open(f'{self.outpath}../failed_runs.log', 'a') as f:
                    f.write(f'{self.outfile}\n')

            # back to working directory
            os.chdir(self.cwd)

            # remove temporary directory
            self.tmpdirobj.cleanup()


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# IMAGING
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i','--inpath',metavar='Absolute path of input data',type=str,default=None, help='Absolute path of input data')
parser.add_argument('-o','--outpath',metavar='Absolute path to save outputs',type=str,default=None, help='Absolute path to save outputs')
parser.add_argument('--hdf5in',metavar='Input hdf5 file',type=str,default=None, help='Input hdf5 file name')
parser.add_argument('--hdf5out',metavar='Output hdf5 file',type=str,default=None, help='Output hdf5 file name')
parser.add_argument('--offset',metavar='ID number offset',type=int,default=0, help='ID number offset')
parser.add_argument('--nproc',metavar='Number of processors to use',type=int,default=0, help='Number of processors to use')
parser.add_argument('--uvfits', help='Input UVFITS file name')
parser.add_argument('--param_search_id', default=[], help='parameter id to search')

inpath = parser.parse_args().inpath
outpath = parser.parse_args().outpath
hdf5in = parser.parse_args().hdf5in
hdf5out = parser.parse_args().hdf5out
offset = parser.parse_args().offset
nproc = parser.parse_args().nproc
uvfits = parser.parse_args().uvfits
param_search_id = parser.parse_args().param_search_id

param_range = []
if param_search_id != 'all':
    param_search_id = param_search_id.split(',')
    for param in param_search_id:
        param_range.append(int(param))

if inpath[-1] != '/':
    inpath += '/'
if outpath[-1] != '/':
    outpath += '/'

if nproc < 1:
    nproc = multiprocessing.cpu_count()
    print('N processor: %s'%(nproc))

hideprint = False

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Convert input parameter table (hdf5) to dictionary
def hdf5_to_pset(hdf5, idrange=[]):
    '''
    idrange: Specify ids range, empty if all (use dataframe INDEX, not id!!)
    '''
    # Load DataFrame, parameters defined above will be override
    print(f'\n... LOADING DATAFRAME {os.path.basename(hdf5)} ...')
    psets = pd.read_hdf(hdf5, 'parameters')
    psets['uvfits'] = [uvfits]*len(psets)
    psets['lmtcal'], psets['jcmtcal'] = False, False
    psets.loc[psets.ref_type == 'False', 'ref_type'] = False

    if len(idrange) > 0:
        psets = psets.take(idrange)

    # Convert to dictionary
    psets = psets.to_dict(orient='records')
    return psets


def main(hdf5in, hdf5out, idrange=[], offset=0):

    # Create outdir if doesn't exist
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    if not os.path.exists('./hdf5s'):
        os.makedirs('./hdf5s')

    # Build DataFrame
    parms = hdf5_to_pset(hdf5in, idrange=idrange)
    psets = paramsurvey.params.product(parms, infer_category=False)

    print(f'\n... IDENTIFIED {len(psets)} PARAMETER COMBINATIONS ...')

    # Remove 'normalized' string from uvfits not normalized
    mask = psets['normalized'] == False
    psets.loc[mask, 'uvfits'] = psets.loc[mask, 'uvfits'].replace({'normalized_10s': '10s'}, regex=True)

    # Drop duplicated rows
    len0 = len(psets)
    psets.loc[psets.deblurr == False, ['ref_type', 'ref_scale']] = 'False', 1.
    psets.loc[psets.psd_noise == False, ['psd_a', 'psd_u0', 'psd_b', 'psd_c']] = 1., 1., 1., 1.
    psets = psets.drop_duplicates(ignore_index=True)
    len1 = len(psets)

    # Add epoch to psets
    psets['epoch'] = np.array([psets['uvfits'][i].split('_')[1] for i in range(len(psets))])

    # Add ids to psets
    ncomb = len(psets) // len(set(psets['uvfits']))
    print(ncomb)

    print(f'\n... {len0-len1} COMBINATIONS WERE DUPLICATED, FINAL NUMBER IS {len(psets)} ...')

    # Save Pandas DataFrame
    psets.to_hdf(hdf5out, 'parameters', mode='w', complevel=9, format='table')

    psets = hdf5_to_pset(hdf5out, idrange=[])


####################
#   Parallelization
####################
    with multiprocessing.Pool(processes=nproc) as pool:
        pool.map(_run_survey, psets)

###################
# Helper functions
###################

def _run_survey(pset):
    with HiddenPrints(hideprint=hideprint):
        survey = FitCleanSurvey(pset, inpath, outpath)
        survey.run()

######
# RUN
######

if __name__ == '__main__':
    start = time.time()
    main(hdf5in, hdf5out, idrange=param_range, offset=offset)
    #main(parms, offset=offset)
    print('\n That took {} minutes'.format((time.time() - start)/60.))
