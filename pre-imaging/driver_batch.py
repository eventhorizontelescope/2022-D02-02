#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This is an example driver of the standardized pre-imaging calibration function.
"""
import preimcal
import ehtim as eh
import argparse
import os

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', "--infile",  metavar='input uvfits file', type=str, required=True, nargs="+",
                    help='input uvfits file name ')

parser.add_argument('--overwrite', action='store_true', default=False,
                    help='If specified, the input file will be overwritten.')

parser.add_argument('--is_normalized', action='store_true', default=False,
                    help='If the input uvfits is flux-normalized, please DO NOT FORGET TO SPECIFY THIS OPTION.')

parser.add_argument('--is_deblurred', action='store_true', default=False,
                    help='If the input uvfits is deblurred, please DO NOT FORGET TO SPECIFY THIS OPTION.')

parser.add_argument("--nproc", metavar="nproc", type=int, default=1,
                    help="Number of processors. Default=1.")

parser.add_argument('--lmtcal', action='store_true', default=False,
                    help='If specified, do LMT cal [NOT RECOMMENDED TO DO]')

parser.add_argument('--lmtfwhm',  metavar='LMTcal FWHM', type=float, default=60.,
                    help='FWHM size in uas for the LMT cal.')

parser.add_argument('--jcmtcal', action='store_true', default=False,
                    help='If specified, do JCMT cal [NOT RECOMMENDED TO DO.]')

parser.add_argument("--tint", metavar='Integ Time', type=float, default=-1,
                    help='Integration Time in seconds. Negative value means no averaging, 0 means scan-averaging.')

parser.add_argument("--syserr", metavar='frac syserr', type=float, default=-1,
                    help='Fractional systematic errors to be added to visibilities in quadrature.')

parser.add_argument('--reftype', metavar='ref type', type=str, default="quarter1",
                    help='Refractive noise table (quarter1, quarter2, dime). If else specified, do not add any budgets. Default to quarter1.')

parser.add_argument('--refscale', metavar='ref scale', type=float, default=1.0,
                    help='Scaling factor of the refractive noise. Default=1.0.')

parser.add_argument('--deblurr', action='store_true', default=False,
                    help='If specified, the output will be deblurred.')

parser.add_argument('--psd_noise', action='store_true', default=False,
                    help='If specified, the PSD noise will be added.')

parser.add_argument('-a', metavar='PSD a', type=float, default=0.018,
                    help='PSD Noise Param: Scaling factor in Jy. Typically (0, 0.1) Jy.')

parser.add_argument('-u0', metavar='PSD u0', type=float, default=2.0,
                    help='PSD Noise Param: Break location of the broken power law in Glambda. Typically (0, 2).')

parser.add_argument('-b', metavar='PSD b', type=float, default=2.5,
                    help='PSD Noise Param: power-law index at longer baselines than u0. Typically (0, 6).')

parser.add_argument('-c', metavar='PSD c', type=float, default=2.0,
                    help='PSD Noise Param: power-law index at shorter baselines than u0. Typically (1.5, 2.5).')

args = parser.parse_args()


if args.reftype.lower() not in ["quarter1", "quarter2", "dime"]:
    reftype = None
else:
    reftype = args.reftype.lower()

for infile in args.infile:
    print("*** Processing %s ***" % (infile))
    if args.overwrite:
        outfile = infile
    else:
        infileh, infilee = os.path.splitext(infile)
        outfile = infileh+"_precal"+infilee

    obs = eh.obsdata.load_uvfits(infile)
    obs_cal = preimcal.preim_pipeline(
        # input observational data
        inputobs=obs,
        # is the input normalized? THIS IS A REALLY IMPORTANT ARGUMENT SO PLEASE BE CAREFUL
        is_normalized=args.is_normalized,
        # is the input deblurred? THIS IS A REALLY IMPORTANT ARGUMENT SO PLEASE BE CAREFUL
        is_deblurred=args.is_deblurred,
        # number of process
        nproc=args.nproc,
        # if needed you can run the LMT cal though unrecommended,
        # since we already distuributed LMT-caled data.
        do_LMTcal=args.lmtcal,
        LMTcal_fwhm=args.lmtfwhm,  # FWHM in uas
        # if needed you can run the JCMT cal though unrecommended,
        # since we already distuributed JCMT-caled data.
        do_JCMTcal=args.jcmtcal,
        # integration time in seconds
        tint=args.tint,
        # fractional systematic errors
        syserr=args.syserr,
        # refractive noise budgets.
        #   The scattering subteam's recommendation is "quater1".
        ref_optype=reftype,
        ref_scale=args.refscale,  # the scaling factor
        # deblurr data or not
        do_deblurr=args.deblurr,
        # add a broken power-law budget for the intra-day variations
        do_psd_noise=args.psd_noise,
        a=args.a,
        u0=args.u0,
        b=args.b,
        c=args.c
    )
    obs_cal.save_uvfits(outfile)
