#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Example driver for the SMILI Sgr A* Stokes I Imaging Pipeline for EHT observations in April 2017

Authors: The Event Horizon Telescope Collaboration et al.
Date: May 10, 2022
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2022c, ApJL, 930, L14 (Sgr A* Paper III)
Data Product Code: 2022-D02-01

Brief Description:
The script is an example driver of the SMILI Sgr A* Stokes I Imaging Pipeline
(smili_pipeline.py) for EHT observations in April 2017 attached in the
same directory. For more detail instructions for smili_pipeline.py
please read the help document associated in the imaging script
"python smili_pipeline.py --help".

This sript reconstructs images on a represent parameters from calibratd uvfits files released in the Data Product Release
(2022-D02-01; see also Sgr A* Paper II). You can run it by
"python example_driver.py --uvfitsdir <xxxx/uvfits> --nproc <number of the processors>."
For details, please have a look at the help document by
"python example_driver.py --help".

References:
 - EHT Collaboration Data Portal Website:
   https://eventhorizontelescope.org/for-astronomers/data
 - The Event Horizon Telescope Collaboration, et al. 2022a, ApJL, 930, L12 (Sgr A* Paper I)
 - The Event Horizon Telescope Collaboration, et al. 2022b, ApJL, 930, L13 (Sgr A* Paper II)
 - The Event Horizon Telescope Collaboration, et al. 2022c, ApJL, 930, L14 (Sgr A* Paper III)
 - The Event Horizon Telescope Collaboration, et al. 2022d, ApJL, 930, L15 (Sgr A* Paper IV)
 - The Event Horizon Telescope Collaboration, et al. 2022e, ApJL, 930, L16 (Sgr A* Paper V)
 - The Event Horizon Telescope Collaboration, et al. 2022f, ApJL, 930, L17 (Sgr A* Paper VI)
 - Akiyama, K., Ikeda, S., Pleau, M., et al. 2017a, ApJ, 838, 1
 - Akiyama, K., Kuramochi, K., Ikeda, S., et al. 2017b, AJ, 153, 159
 - Moriyama, K., Akiyama, K., Cho, I., et al. 2019, Zenodo (SMILI version 0.2.0)
 - SMILI: https://github.com/astrosmili/smili
"""

# -------------------------------------------------------------------------------
# Information of Authors
# -------------------------------------------------------------------------------
__author__ = "The Event Horizon Telescope Collaboration et al."
__copyright__ = "Copyright 2022, the Event Horizon Telescope Collaboration et al."
__license__ = "GPL version 3"
__version__ = "1.0"
__date__ = "May 10 2022"


# -------------------------------------------------------------------------------
# Modules
# -------------------------------------------------------------------------------
import tqdm
import os
import argparse
import ray
import pandas as pd
import pickle
from glob import glob
pickle.HIGHEST_PROTOCOL = 4


# -------------------------------------------------------------------------------
# Loading command line arguments
# -------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--uvfitsdir', metavar='Input uvfits directory', type=str, default=None, required=True,
                    help='Input uvfits directory of the data product release (e.g., 2022-XX-XXXX)')
parser.add_argument('--parameter', metavar='Input parameter table', type=str, default=None, required=True,
                    help='Input parameter table')
parser.add_argument('-o', '--outputdir', metavar='Output fits file directory', type=str, default='./smili_reconstructions',
                    help='Output directory. The default is "./smili_reconstructions".')
parser.add_argument('--pipeline', metavar='Pipeline script', type=str, default='smili_pipeline.py',
                    help='The filename of the pipeline script. The default is "smili_pipeline.py".')
parser.add_argument('--nproc', metavar='Number of parallel processors. Default: all cpus.', type=int, default=None)

# -------------------------------------------------------------------------------
# Pipeline Class
# -------------------------------------------------------------------------------
class PipelineDriver(object):
    params_default = dict(
        id=0,
        nx=80,
        pol="PI",
        normalized=True,
        tint=60., wint=-1,
        syserr=0.02,
        l1_lambda=1.0, l1_noise=0.1,
        tv_lambda=1e+3, tsv_lambda=1e+3,
        gs_lambda=-1,
        timeframe="static", rt_lambda=-1, ri_lambda=-1,
        wamp=1.,
        shift_type="com1",
        lc_type="zbl", extflux=0.,
        init_im="tophat", init_size=60.,
        prior_im="tophat", prior_size=100.,
        snrcutoff=0.,
        lmtcal=False, lmtfwhm=-1,
        deblurr=True, ref_type="quarter1", ref_scale=1.0,
        psd_noise=True, psd_a=0.02, psd_u0=2.0, psd_b=2.0, psd_c=2.0
    )

    def __init__(
        self, inuvfits, paramtable,
        filehead="dataset", outdir="./output",
        pipeline="smili_pipeline.py",
    ):
        self.paramtable = paramtable
        self.Npar = len(paramtable.id)
        self.outdir = outdir
        self.filehead = filehead
        self.pipeline = pipeline

        # create output directory
        if not os.path.isdir(outdir):
            os.makedirs(outdir, exist_ok=True)

    def get_command(self, idx):
        # sanity check
        if idx >= self.Npar or idx < 0:
            raise ValueError(
                "The input index is not in the appropriate range")

        # input parameter
        params_input = self.paramtable.loc[idx, :].to_dict()

        # update pipeline parameters
        params = self.params_default.copy()
        params_keys = params.keys()
        for key in params_input.keys():
            if key in params_keys:
                params[key] = params_input[key]

        # get the output file name
        if params["deblurr"]:
            sctcode = "dsct"
        else:
            sctcode = "sct"
        outfits = os.path.join(
            self.outdir,
            "{head}_{id:08d}_{sctcode}.fits".format(
                head=self.filehead,
                id=params["id"],
                sctcode=sctcode
            )
        )
        # generate command
        command = []
        command.append("python")
        command.append(self.pipeline)
        command.append("-i %s" % (" ".join(uvfitsfiles)))
        command.append("-o %s" % (outfits))

        command.append("--nproc 1")

        # Image FoV
        command.append("--nx %d" % (params["nx"]))
        command.append("--pol %s" % (params["pol"]))

        # Flux Normalization
        if params["normalized"]:
            command.append("--is_normalized")

        # Time averaging & Reweighting
        command.append("--tint %g" % (params["tint"]))
        command.append("--wint %g" % (params["wint"]))

        # Systematic Error
        command.append("--sys %g" % (params["syserr"]))

        # Initital Image
        command.append("--init_im %s" % (params["init_im"]))
        if params["init_im"] in ["tophat", "gaussian"]:
            command.append("--init_im_size %g" % (params["init_size"]))

        # Prior Image
        command.append("--prior_im %s" % (params["prior_im"]))
        if params["prior_im"] in ["tophat", "gaussian"]:
            command.append("--prior_size %g" % (params["prior_size"]))

        # Static Image Regularization
        command.append("--lambl1 %g" % (params["l1_lambda"]))
        command.append("--l1_noise %g" % (params["l1_noise"]))
        command.append("--lambtv %g" % (params["tv_lambda"]))
        command.append("--lambtsv %g" % (params["tsv_lambda"]))
        command.append("--lambgs %g" % (params["gs_lambda"]))

        # Dynamical Regularization
        if (params["timeframe"] == "static") or (params["timeframe"] is None) or (params["timeframe"] is False):
            command.append("--lambrt -1")
            command.append("--lambri -1")
        else:
            command.append("--time_frame %s" % (params["timeframe"]))
            command.append("--lambrt %g" % (params["rt_lambda"]))
            command.append("--lambri %g" % (params["ri_lambda"]))

        # Data Weight
        command.append("--wamp %g" % (params["wamp"]))

        # SNR cutoff
        command.append("--snrcutoff %g" % (params["snrcutoff"]))

        # Shift
        command.append("--shifttype %s" % (params["shift_type"]))

        # Total Flux
        command.append("--lc_type %s" % (params["lc_type"]))
        command.append("--extflux %g" % (params["extflux"]))

        # LMT Calibration
        if params["lmtcal"] and params["lmtfwhm"] > 0:
            command.append("--lmtfwhm %g" % (params["lmtfwhm"]))
        else:
            if params["lmtcal"]:
                raise ValueError("lmtcal and lmtfwhm is not consistent.")
            command.append("--lmtfwhm -1")

        # Scattering
        if params["deblurr"]:
            command.append("--deblurr")

        if params["ref_type"].lower() != "none":
            if params["ref_scale"] <= 0:
                raise ValueError("ref_scale value is negative or zero.")
            command.append("--reftype %s" % (params["ref_type"]))
            command.append("--refscale %g" % (params["ref_scale"]))
        else:
            command.append("--reftype none")
            command.append("--refscale 1.0")

        # PSD Noise
        if params["psd_noise"]:
            if (params["psd_a"] < 0) or (params["psd_u0"] < 0) or (params["psd_b"] < 0) or (params["psd_c"] < 0):
                raise ValueError(
                    "invalid psd values are specified despite psd_noise=True.")
            command.append("--a %g" % (params["psd_a"]))
            command.append("--u0 %g" % (params["psd_u0"]))
            command.append("--b %g" % (params["psd_b"]))
            command.append("--c %g" % (params["psd_c"]))
        else:
            command.append("--a -1")
            command.append("--u0 -1")
            command.append("--b -1")
            command.append("--c -1")


        # concat all arguments
        command = " ".join(command)
        return command


if __name__ == "__main__":
    args = parser.parse_args()

    # include uvfits files (the convention of file name is based on ER6 data set)
    uvfitsfile_list  = sorted(glob(args.uvfitsdir+"/ER6_SGRA_2017_*_lo_*_netcal-LMTcal-norm_StokesI.uvfits"))
    commands = []
    for uvfitsfile in uvfitsfile_list:
        uvfitsfiles = [uvfitsfile, uvfitsfile.replace("_lo_", "_hi_")]
        epoch = uvfitsfile.split("/")[-1].split("_")[3]
        if "hops" in uvfitsfile:
            # SGRA hops
            if "SGRA" in uvfitsfile:
                label="hops_%s_LO+HI"%(epoch)
            # synthetic data
            else:
                label = uvfitsfile.split("/")[-1].split("hops_%s_"%(epoch))[1].split("_LO")[0]+"_%s_LO+HI"%(epoch)
        # SGRA casa
        if "casa" in uvfitsfile:
            label="casa_%s_LO+HI"%(epoch)

        paramtable = pd.read_csv(args.parameter)

        outdir = args.outputdir + "/" + label
        filehead = label
        pipeline = args.pipeline
        nproc = args.nproc

        # initialize the pipeline driver
        pipe = PipelineDriver(
            inuvfits=uvfitsfiles,
            paramtable=paramtable.reset_index(),
            filehead=filehead,
            outdir=outdir,
            pipeline=pipeline,
        )
        Npar = len(paramtable.id)
        commands += [pipe.get_command(idx)
                    for idx in tqdm.tqdm(range(Npar), desc="Generate Commands")]
        del paramtable
        del pipe

#-------------------------------------------------------------------------------
# Running the pipeline
#-------------------------------------------------------------------------------

#commands = list(set(commands))
#for command in commands:
#    os.system(command)
    ray.init(num_cpus=nproc)

    @ray.remote
    def run_command(command):
        #print(command)
        return os.system(command)

    commands = list(set(commands))
    output = [run_command.remote(command)
              for command in tqdm.tqdm(commands, desc="Ray Remote")]

    ray.get(output)
