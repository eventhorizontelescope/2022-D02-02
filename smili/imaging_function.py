#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ehtim as eh
import numpy as np
import pandas as pd
import pickle
import scipy.interpolate as interp
import ehtim.imaging.dynamical_imaging as di
import ehtim.scattering.stochastic_optics as so
from smili import uvdata, imdata, util
import copy
import astropy.time as at
import ehtplot
import random
import os

################################################################################
# Prior functions
################################################################################
def set_l1prior(prior_image, prior_size, dx_uas=2, nx=80, reset_xyref=True):
    '''
    Make prior image with arguments in a imaging pipeline.
    Args:
        prior_image: Prior image structure (tophat, gaussian, or input fits file name)
        prior_size: Prior size (uas)
    '''
    if ".fits" in prior_image:
        l1_prior = imdata.IMFITS(dx=dx_uas, nx=nx, angunit="uas")
        if reset_xyref:
            l1_prior.header["nxref"] =int(l1_prior.header["nxref"])
            l1_prior.header["nyref"] =int(l1_prior.header["nyref"])

        l1_prior = l1_prior.cpimage(imdata.IMFITS(prior_image, "ehtim"))
    elif prior_image == "gaussian":
        l1_prior = gen_prior_gauss(
            fwhm=prior_size, epsmin=-1, dx=dx_uas, nx=nx, reset_xyref=reset_xyref).copy()
    elif prior_image == "tophat":
        l1_prior = gen_prior_tophat(
            radius=prior_size/2., epsmin=-1, dx=dx_uas, nx=nx,gauss_fwhm=20., reset_xyref=reset_xyref).copy()
    else:
        raise ValueError("Information of prior_image is incorrect.")
    return l1_prior

def gen_prior_gauss(x0=0., y0=0., fwhm=40., epsmin=-1, dx=2, nx=80, reset_xyref=True):
    '''
    This function provides a Gaussain image with a specified FWHM in uas.
    '''
    # create a blank image
    outimage = imdata.IMFITS(dx=dx, nx=nx, angunit="uas")
    if reset_xyref:
        outimage.header["nxref"] =int(outimage.header["nxref"])
        outimage.header["nyref"] =int(outimage.header["nyref"])

    outimage = outimage.add_gauss(majsize=fwhm, overwrite=True)
    if epsmin > 0:
        outimage.data += outimage.peak()*epsmin
        outimage.data /= outimage.totalflux()
    return outimage

def gen_prior_tophat(
        radius=40.,
        totalflux=1.,
        epsmin=-1,
        dx=2, nx=80,
        gauss_fwhm=10,
        reset_xyref=True):
    '''
    This function provides a top hat image with a specified radius in uas.
    '''
    image = imdata.IMFITS(dx=dx, nx=nx, angunit="uas")
    if reset_xyref:
        image.header["nxref"] =int(image.header["nxref"])
        image.header["nyref"] =int(image.header["nyref"])
    x, y = image.get_xygrid(twodim=True, angunit="uas")
    rr = np.sqrt(x**2+y**2)
    idx = rr <= radius
    jnu = copy.deepcopy(rr)
    jnu[:] = 0
    jnu[idx] = 1
    image.data[0, 0, :, :] = jnu
    image.data = image.data*totalflux/image.totalflux()
    if epsmin > 0:
        image.data += image.peak()*epsmin
        image.data /= image.totalflux()
    image = image.convolve_gauss(gauss_fwhm)

    return image

################################################################################
# Align functions in imaging iteration
################################################################################
def edit_movie(movie, imprm, cbeamprm, shifttype=None):

    '''
    Movie editing (shifting, blurring)
    Args:
        movie: imdata.MOVIE object
        imprm: imaging arguments
        cbeamprm: beam parameters
        shifttype: type of movie centering ("nxcorr","peak","hough")
    Return:
        imdata.MOVIE object
    '''

    movie_edited = movie.copy()
    if "vistable" not in list(imprm.keys()):
        lc_frm = movie.make_lightcurve()["flux"]

        for it in range(movie_edited.Nt):
            movie_edited.images[it] = movie_edited.images[it].convolve_gauss(
                scale=1./4, **cbeamprm).copy()

        # Shifting images such that the ring or dominant structure comes to the center
        alpha = np.int32(shifttype.lower().split("com")[1])
        movie_edited.images = align_images_com(
            movie_edited.images, save_totalflux=False, alpha=alpha)[0]

        # Adjust light curve
        for it in range(movie_edited.Nt):
            movie_edited.images[it].data *= lc_frm[it] / \
                movie_edited.images[it].totalflux()

    return movie_edited

def align_images_com(images, alpha, save_totalflux=False):
    '''
    Align frame images using com centering
    Args:
        images: list of imdata.IMFITS object
        save_totalflux: if specified, totalflux is conserved
        comshift: if specified, the image center is set to
                  that of com origin of the mean image
    Return:
        outimages: list of aligned imdata.IMFITS object
        meanimage: mean image of outimages
    '''
    # number of images
    Nimage = len(images)

    # Shift Images
    totalfluxes = [image.totalflux() for image in images]
    outimages = [image.peakshift() for image in images]
    meanimage = copy.deepcopy(outimages[0])
    meanimage.data[:] = np.mean(np.asarray(
        [outimages[i].data[:] for i in range(Nimage)]), axis=0)

    # get ceter of the mass
    compos = meanimage.compos(alpha=alpha)
    meanimage = meanimage.refshift(**compos)

    # com shift
    outimages = [image.refshift(**compos) for image in outimages]

    if save_totalflux:
        meanimage.data *= np.sum(totalfluxes)/Nimage/meanimage.totalflux()
        for i in range(Nimage):
            outimages[i].data *= totalfluxes[i]/outimages[i].totalflux()

    return outimages, meanimage

################################################################################
# Initial condition of imagings
################################################################################
def set_initialmovie_pipeline(init_im, init_im_size, time_frame, uvfits_file, dx_uas, nx, timetype="uniform"):
    '''
    Set initial movie from each argument in imaging pipeline
    Args:
        init_im: initial image structure (tophat, gaussian, or input fits file name)
        init_im_size: initial image size (uas)
        time_frame: information for a time frame of a movie (time interval (minutes), or time table name
        uvfits_file: uvfits file name
    Return:
        imdata.MOVIE object
    '''

    init_image = set_initialimage(
        init_im, init_im_size, uvfits_file, dx_uas, nx)
    initmovie = set_initialmovie(init_image, time_frame, uvfits_file, timetype=timetype)
    return initmovie


def set_initialimage(init_im_key, init_im_size, uvfits_file, dx_uas=2, nx=80):
    '''
    Set intial image
    Args:
        init_im_key: initial image structure (tophat, gaussian, or input fits file name)
        init_im_size: size of a tophat or gaussian
        uvfits_file: uvfits file name
    Return:
        imdata.IMFITS object
    '''

    if ".fits" in init_im_key:
        image = imdata.IMFITS(dx=dx_uas, nx=nx, angunit="uas").cpimage(
            imdata.IMFITS(init_im_key, uvfits=uvfits_file, imfitstype="ehtim"))
    elif init_im_key == "tophat" or init_im_key == "gaussian":
        image = imdata.IMFITS(
            dx=dx_uas, nx=nx, angunit="uas", uvfits=uvfits_file)
        init_im_size = init_im_size
        if init_im_key == "gaussian":
            image = image.add_gauss(
                totalflux=1., majsize=init_im_size, overwrite=True)
        elif init_im_key == "tophat":
            image = gen_prior_tophat(
                radius=init_im_size/2., epsmin=-1, dx=dx_uas, nx=nx)
        del init_im_size

        image.header["nxref"] = int(image.header["nxref"])
        image.header["nyref"] = int(image.header["nyref"])
    return image


def set_initialmovie(image, time_frame_key, uvfitsfile, timetype="uniform"):
    '''
    Make initial movie whose frame images are the same as input image
    Args:
        image: imdata.IMFITS object
        time_frame_key: information for a time frame of a movie (time interval (minutes), or time table name)
        uvfitsfile: file name of uvfits
        timetype: structure of the movie time-table. if "uniform" is uniform time frame,
                       while "non-uniform" does not includes time frame without visibility data
    Return:
        imdata.MOVIE object
    '''

    # Set dateobs header for ehtim inclusion
    obs = eh.obsdata.load_uvfits(uvfitsfile)
    dateobs = at.Time(obs.mjd, format="mjd")
    dateobs = "%04d-%02d-%02d" % (dateobs.datetime.year,
                                  dateobs.datetime.month, dateobs.datetime.day)

    angunit = image.angunit
    image.header["object"] = obs.source
    image.header["x"] = obs.ra * 12 * util.angconv(angunit, "deg")
    image.header["y"] = obs.dec * util.angconv(angunit, "deg")
    image.header["f"] = obs.rf
    image.header["dateobs"] = dateobs

    # Set time frame
    if ".csv" in time_frame_key:
        initmovie  = imdata.MOVIE(timetable=time_frame_key)

    else:
        uvfits = uvdata.UVFITS(uvfitsfile).select_stokes("PI")
        vtable = uvfits.make_vistable()
        utcst  = vtable.get_utc()[0]
        utced = vtable.get_utc()[-1]

        if time_frame_key == "static":
            dt = at.TimeDelta(utced-utcst)
            timebin = at.Time(np.array([utcst+dt*0.5]))
        else:
            dt = at.TimeDelta(float(time_frame_key)*60, format="sec")
            timebin = at.Time(np.arange(utcst, utced+0.5*dt, dt))

        if timetype=="non-uniform":
            idx_c = []
            tobs = vtable.get_utc().cxcsec
            tmov = timebin.cxcsec
            for it in range(len(tmov)):
                tinit = tmov[it]-dt.sec/2
                tend = tmov[it]+dt.sec/2
                idx = tobs >tinit
                idx &= tobs<tend
                num = len(tobs[idx])
                if num==0:
                    idx_c.append(False)
                else:
                    idx_c.append(True)
            timebin = timebin[idx_c]

        initmovie = imdata.MOVIE(tcen=timebin, tint=dt.sec)

    # Set initial image
    initmovie = initmovie.set_image(image)
    return initmovie


################################################################################
# Light curve functions
################################################################################
def set_lightcurve(lc_type, exflux, uvf_list, baseline_list=[["AA", "AP"],["JC", "SM"]], kind="nearest"):
    '''
    Make lightcurve based on light curve information.
    Args:
        lc_type: Light curve information (zbl, uniform total flux value, or table name)
        exflux: Extened flux (Jy). Light curve (lc) of a movie = original lc -external flux
        uvf_list: list of uvfits name or uvdata.UVFITS
        baseline_list: base line of zero base line (only for lc_type="zbl")
        kind: Interpolation type of light curve
    Return:
        lightcurve object
    '''

    if not isinstance(uvf_list[0], uvdata.UVFITS) and ".uvf" in uvf_list[0].lower():
        uvfits_list = [uvdata.UVFITS(inputuvfile) for inputuvfile in uvf_list]
    else:
        uvfits_list = uvf_list

    if ".csv" in lc_type:
        lc = lightcurve.read_csv(lc_type)

    elif lc_type == "zbl":
        flux_list = []
        i = 0
        for uvfits in uvfits_list:
            vtable_tmp = uvfits.select_stokes("PI").make_vistable()
            if i == 0:
                # Set time reference
                lc_ref = vtable_tmp.make_lightcurve(
                    baseline_list=baseline_list).reset_index(drop=True)
                flux_list.append(lc_ref["flux"])
            else:
                lc_zbl = lc_ref.cplightcurve(vtable_tmp.make_lightcurve(
                    baseline_list=baseline_list), kind=kind).reset_index(drop=True)
                flux_list.append(lc_zbl["flux"])
            i += 1
        lc = copy.deepcopy(lc_ref)
        lc["flux"] = np.mean(flux_list, axis=0)

    else:
        vtable_tmp = uvfits_list[0].select_stokes("PI").make_vistable()
        if baseline_list is not None:
            lc = vtable_tmp.make_lightcurve(baseline_list=baseline_list).reset_index(drop=True)
        else:
            lc = vtable_tmp.make_lightcurve().reset_index(drop=True)
        lc["flux"][:] = np.float64(lc_type)

    lc["flux"] -= exflux
    return lc


def set_lightcurve_movie(movie, lc_type, exflux, uvf_list, baseline_list=[["AA", "AP"], ["JC", "SM"]], kind="nearest"):
    '''
    Make light curve using time-frame information of a movie.

    Args:
        movie: imdata.MOVIE object
        lc_type: Light curve information (zbl, uniform total flux value, or table name)
        exflux: Extened flux (Jy). Light curve (lc) of a movie = original lc -external flux
        uvf_list: list of uvfits name or uvdata.UVFITS
        baseline_list: base line of zero base line (only for lc_type="zbl")
        kind: Interpolation type of light curve
    Return:
        lightcurve object
    '''

    lctable = set_lightcurve(lc_type, exflux, uvf_list, baseline_list, kind)
    #lcarr = make_zbl_lcarr(obs)
    #lcinterp = interp.interp1d(lcarr["time"].flatten(), lcarr["amp"].flatten(),
    #                           fill_value="extrapolate")
    #lcflux = lcinterp(obs.data['time'])
    lc_frm = movie.make_lightcurve()
    lctable = lc_frm.cplightcurve(lctable, kind=kind)
    return lctable
