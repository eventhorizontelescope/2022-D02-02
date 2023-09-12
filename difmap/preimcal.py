#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ehtim as eh
import numpy as np
import pandas as pd
import pickle
import scipy.interpolate as interp
import ehtim.imaging.dynamical_imaging as di
import ehtim.scattering.stochastic_optics as so

################################################################################
# Pre-imaging Calibration
################################################################################


def preim_pipeline(
        inputobs,
        is_normalized=False,
        is_deblurred=False,
        lcarr=None,
        nproc=1,
        do_LMTcal=False,
        LMTcal_fwhm=60.,
        do_JCMTcal=False,
        tint=-1,
        syserr=-1,
        ref_optype="quarter1",
        ref_scale=1.,
        do_deblurr=True,
        do_psd_noise=True,
        a=0.018,
        u0=2.,
        b=2.5,
        c=2.):
    """
    Pre-imaging pipeline function. This function is desined to perform the
    standardized pre-imaging processing for 2017 Stokes Sgr A* data.

    Args:
        inputobs (ehtim.obsdata.Obsdata):
            Input observational data.
        is_normalized (bool, optional):
            If the input data is flux-normalized. Defaults to False.
        is_deblurred (bool, optional):
            If the input data is deblurred. Defaults to False.
        lcarr:
            Light curve table. If not specified, guess from the intra-site
            flux density.
        nproc (int, optional):
            Number of processing for some procedures. Defaults to 1.
        doLMTcal (bool, optional):
            If True, do LMT cal. Not recommended for the parameter survey,
            since it is a bit slow and we already provided LMT-calibrated data.
            Defaults to False.
        LMTcal_fwhm (int, optional):
            FWHM size for the LMT cal. Defaults to 60.
        do_JCMTcal (bool, optional):
            If True, do HCME phase cal. Not recommended for the parameter survey,
            since it is a bit slow and we already provid JCMT calibrated data.
            Defaults to False.
        tint (float, optional):
            Integration time in seconds. Skipped If negative.
            Defaults to -1.
        syserr (float, optional):
            Fractional systematic errors. Skipped If negative.
            Defaults to -1.
        ref_optype (str, optional):
            If not None, add the refractive noise floor. All based on Johnson+18 model.
            You can switch the type of the floor.
            The recommended default model is "quarter1"
                'dime': add a constant value of (10 mJy * scale) in quadrature to all (u,v) data;
                'quarter1': (u,v)-dependent noisefloor values based on the RMS renormalized
                            refractive noise values of a circular Gaussian model;
                'quarter2': (u,v)-dependent noisefloor values based on the averaged RMS renormalized
                            refractive noise values of several synthetic models (ring, crescent, GRMHD, Gaussian).
        ref_scale (float, optional):
            The scaling factor for the refractive noise floor. Defaults to 1.
        do_deblurr (bool, optional):
            If True, data will be deblurred. Defaults to True.
        do_psd_noise (bool, optional):
            If True, the error budgets from the intra-day variability
            will be added based on a broken power-law model, whose parameters
            are given by a, u0, b and c. Defaults to True.
        a (float, optional):
            The scaling factor of the power spectrum amplitude.
            Typically (0, 0.1) Jy.
        u0 (float, optional):
            The location of the break in the broken-power law spectrum,
            in the unit of GLambda. Typically (0, 2) GLambda.
        b (float, optional):
            The power law index at longer baselines than u0. Typically (0,6).
        c (float, optional):
            The power law index at shorter baselines than u0. Typically (1.5, 2.5).

    Returns:
        ehtim.obsdata.Obsdata object: calibrated observational data
    """
    # Copy the input obsdata to the one for the output.
    obs = inputobs.copy()

    # Check the input light curve. If it is not specified, guess it from
    # the intra-site VLBI flux density
    if lcarr is None and (not is_normalized or do_LMTcal):
        print("Pre-imaging processing: Obtaining the light curve from the intra-site baselines")
        lcarr = make_zbl_lcarr(obs)

    # Network calibration if needed
    #   Networkcal should use the lightcurve for large-scale (arcsec-scale)
    #   emission (=lc_netcal) anticipated for intra-site baselines
    #   obs = netcal(uvfile, lcarr, nproc=nproc, smooth=0.9, num_repeat=3)

    # Normalize flux
    #   If the input obsfile is not flux normalized, normalize all visibility
    #   amplitudes with the light curve.
    if not is_normalized:
        print("Pre-imaging processing: Normalizing the total flux density")
        lcinterp = interp.interp1d(lcarr["time"].flatten(), lcarr["amp"].flatten(),
                                   fill_value="extrapolate")
        lcflux = lcinterp(obs.data['time'])
        obs_rs = obs.switch_polrep('circ')
        for field in ['rrvis', 'llvis', 'rlvis', 'lrvis', 'rrsigma', 'llsigma', 'rlsigma', 'lrsigma']:
            obs_rs.data[field] /= lcflux

        obs = obs_rs.copy()

    # Blur data
    if is_deblurred:
        print("Pre-imaging processing: reblurring data.")
        sm = so.ScatteringModel()
        obs = Blur_obs(sm, obs.switch_polrep("stokes"))

    # LMT calibration if specified.
    if do_LMTcal and ("LM" in obs.data["t1"] or "LM" in obs.data["t2"]):
        print("Pre-imaging processing: LMT calibration")
        #   note: since data is flux-nomarlized, this needs static LMT cal.
        obs = LMT_calibration(
            obs,
            lcarr=np.ones_like(lcarr["amp"]),
            caltype="const",
            LMTcal_fwhm=LMTcal_fwhm,
            nproc=nproc
        )

    #  JCMT phase calibration
    if do_JCMTcal and ("JC" in obs.data["t1"] or "JC" in obs.data["t2"]):
        print("Pre-imaging processing: JCMT calibration")
        obs = eh.netcal(
            obs.switch_polrep("circ"),
            1.0, sites=['JC'], method='phase', pol='RRLL', processes=nproc
        )

    # Time averaging
    if tint >= 0:
        print("Pre-imaging processing: time-averaging.")
        obs = obs.avg_coherent(inttime=tint)

    # Add Systematic Errors
    if syserr > 0:
        print("Pre-imaging processing: adding the fractional systematic errors.")
        obs = obs.add_fractional_noise(syserr)

    # add the error of the refractive noise to scattered data
    if (ref_optype is not None) and (ref_optype is not False):
        #   note: this should be added to unblurred but flux-normalized visibilities.
        #         this is why data are flux-normalized but not deblurred
        #         prior to this procedure.
        print("Pre-imaging processing: adding the refractive noise to the flux-normalized data.")
        obs = add_noisefloor_obs(
            obs.switch_polrep("stokes"),
            optype=ref_optype, scale=ref_scale)

    # add noise budgets to scattered data
    add_psd_noise_flag = do_psd_noise and a >= 0. and b > 0 and u0 >= 0 and c >= 0
    if add_psd_noise_flag:
        print("Pre-imaging processing: adding the PSD noise to deblurred data.")
        obs = add_psd_noise(obs.switch_polrep("stokes"), a, u0, b, c)

    # deblurr data if the user needs
    if do_deblurr:
        # deblur data
        print("Pre-imaging processing: deblurring data")
        sm = so.ScatteringModel()
        obs = sm.Deblur_obs(obs.switch_polrep("stokes"))

    # If the input obsdata is not flux normalized, re-scale the total flux.
    if not is_normalized:
        print("Pre-imaging processing: rescaling the total flux density.")
        obs_rs = obs.switch_polrep('circ')
        lcflux = lcinterp(obs_rs.data['time'])
        for field in ['rrvis', 'llvis', 'rlvis', 'lrvis', 'rrsigma', 'llsigma', 'rlsigma', 'lrsigma']:
            obs_rs.data[field] *= lcflux
        obs = obs_rs.copy()

    return obs.switch_polrep("stokes")


################################################################################
# Noise modeling function
################################################################################


def add_psd_noise(obs, a=0.5, u0=1, b=2, c=0, debias=False):
    '''
    Increase errors by specified fractional values of amplitudes

    Args:
        obs: Obsdata object
        a: parameter controls power spectrum amplitude. [typically (0, 0.5) Jy]
        u0 [Glambda]: controls the location of he break in the power law spectrum (1, 5) Glambda
        b: power law index of error at the long baselines (2, 6)
        c: power law index of error at the short baselines (0, 3)
    Returns:
        Obsdata object with the addtional noise
    '''

    # Extract visibility amplitudes
    # Switch to Stokes for graceful handling of circular basis products missing RR or LL
    # amp = obs.switch_polrep('stokes').unpack('amp', debias=debias)['amp']
    uvdist = obs.switch_polrep('stokes').unpack('uvdist')['uvdist']

    out = obs.copy()
    for sigma in ['sigma1', 'sigma2', 'sigma3', 'sigma4']:
        try:
            field = obs.poldict[sigma]
            out.data[field] = (obs.data[field]**2
                               + np.abs(a**2 * ((4./u0)**c / (1. + (4./u0)**(b+c)))**(-1) * (uvdist/(u0*1e+9))**c / (1. + (uvdist/(u0*1e+9))**(b+c))))**0.5
        except KeyError:
            continue

    return out


################################################################################
# Precalibration
################################################################################
def netcal(uvfits, lcarr, nproc=1, smooth=0.5, num_repeat=2):
    '''
    Network calibration
    Args:
        uvfits: eh.obsdata.Obsdata object or uvfits file name
        lcarr: lightcurve array
        nproc: number of parallel processors
        smooth: smoothing parameter
    Return:
        eh.obsdata.Obsdata object
    '''
    if type(uvfits) != eh.obsdata.Obsdata:
        obs = eh.obsdata.load_uvfits(uvfits, polrep='circ')
    else:
        obs = uvfits.copy().switch_polrep('circ')

    spl = make_spline(lcarr, smooth)

    # univariate spline interpolation
    obs_nc = eh.netcal(obs, spl, processes=nproc, gain_tol=1.0, pol='RRLL')
    for repeat in range(num_repeat):
        obs_nc = eh.netcal(obs_nc, spl, processes=nproc,
                           gain_tol=1.0, pol='RRLL')

    return obs_nc


def make_spline(lcarr, smooth=0.9):
    '''
    Make spline object for an input light curve
    Args:
        lcarr: light curve array
        smooth: smoothing parameter
    Return:
        interp.UnivariateSpline object
    '''
    spl = interp.UnivariateSpline(lcarr["time"], lcarr["amp"], ext=3)
    spl.set_smoothing_factor(smooth)
    return spl


def LMT_calibration(obs, lcarr, caltype="const", LMTcal_fwhm=60.0, nproc=1):
    '''
    LMT precalibration
    Args:
        obs: eh.obsdata.Obsdata object
        lcarr: lightcurve object
        caltype: if "const" (zbl) calibration with constant (time variable) gaussian
        LMTcal_fwhm: FWHM of gaussian (uas)
    '''
    fov_uas = 150.0
    nx = 32

    totalflux = np.median(lcarr)
    obs_selfcal = obs.switch_polrep('stokes')
    if caltype == "const":
        #   Set a circular Gaussian with a FWHM of 60 uas as the model image
        gausspriormodel = eh.image.make_square(
            obs, int(nx), fov_uas*eh.RADPERUAS)           # Create empty image
        gausspriormodel = gausspriormodel.add_gauss(
            totalflux, (LMTcal_fwhm*eh.RADPERUAS, LMTcal_fwhm*eh.RADPERUAS, 0, 0, 0))  # Add gaussian
        #   Self-calibrate amplitudes
        for repeat in range(3):
            caltab = eh.self_cal.self_cal(
                obs_selfcal.flag_uvdist(uv_max=2e9),
                gausspriormodel,
                sites=['LM', 'LM'],
                method='vis',
                ttype='nfft',
                processes=nproc,
                caltable=True,
                gain_tol=1.0)
            obs_selfcal = caltab.applycal(
                obs_selfcal, interp='nearest', extrapolate=True)
    elif caltype == "lc":
        # Variable LMT calibration
        spl = make_spline(lcarr)
        obs.add_scans()
        obs_split = obs_selfcal.split_obs()
        for j in range(len(obs_split)):
            zbl = spl(obs_split[j].data['time'][0])
            gausspriormodel = eh.image.make_square(
                obs, int(nx), fov_uas*eh.RADPERUAS)           # Create empty image
            gausspriormodel = gausspriormodel.add_gauss(
                zbl, (LMTcal_fwhm*eh.RADPERUAS, LMTcal_fwhm*eh.RADPERUAS, 0, 0, 0))  # Add gaussian
            for repeat in range(1):
                caltab = eh.self_cal.self_cal(
                    obs_selfcal.flag_uvdist(uv_max=2e9),
                    gausspriormodel,
                    sites=['LM', 'LM'],
                    method='vis',
                    ttype='nfft',
                    processes=nproc,
                    caltable=True,
                    gain_tol=1.0)
                obs_split[j] = caltab.applycal(
                    obs_split[j], interp='nearest', extrapolate=True)
            # Save calibrated uvfits data sets
            obs_selfcal = di.merge_obs(obs_split)
    else:
        obs_selfcal = obs.copy()

    return obs_selfcal


def make_zbl_lcarr(obs):
    '''
    Make light curve array using visibility amplitudes of intrasite baselines
    '''
    AAAP = obs.unpack_bl(site1="AA", site2="AP", fields=["amp"])
    JCSM = obs.unpack_bl(site1="JC", site2="SM", fields=["amp"])
    if len(AAAP) > 1 and len(JCSM) == 0:
        # print("Get the light curve from AA-AP")
        return AAAP
    elif len(AAAP) == 0 and len(JCSM) > 1:
        # print("Get the light curve from JC-SM")
        return JCSM
    elif len(AAAP) > 1 and len(JCSM) > 1:
        # print("Get the light curve from AA-AP and JC-SM")
        # for merge time of spline interpolation
        JCSM["time"][:, 0] += 1e-3
        lcarr = np.concatenate([AAAP, JCSM])
        time = lcarr["time"][:, 0]
        amp = lcarr["amp"][:, 0]
        idx = np.argsort(time)
        lcarr["time"][:, 0] = time[idx]
        lcarr["amp"][:, 0] = amp[idx]
        return lcarr
    else:
        raise ValueError("There is no intrasites of AA-AP and JC-SM")

################################################################################
# Functions for scattering prescription
################################################################################


def add_noisefloor_obs(obs, optype="quarter1", scale=1.0):
    """ Function to add noisefloor to Obsdata.
    Requirements: numpy, scipy, pandas, ehtim;
    'obs' should be an ehtim.obsdata.Obsdata object
    Options for noisefloor include:
    'dime': add a constant value of (10 mJy * scale) in quadrature to all (u,v) data;
    'quarter1': (u,v)-dependent noisefloor values based on the RMS renormalized
                refractive noise values of a circular Gaussian model;
    'quarter2': (u,v)-dependent noisefloor values based on the averaged RMS renormalized
                refractive noise values of several synthetic models (ring, crescent, GRMHD, Gaussian).
    """
    epochs = {57849: '3598', 57850: '3599'}
    epoch = epochs[obs.mjd]
    if optype == "quarter1":
        # pickle file for the noise floor template
        infile = "./obs_scatt_std_gauss_scaled_%s.pickle" % (epoch)
        obs_new = add_noise_floor_obs_quarter(
            obs, infile=infile, scale=scale / 2.27)
    elif optype == "quarter2":
        # pickle file for the noise floor template
        infile = "./obs_scatt_std_avg_scaled_%s.pickle" % (epoch)
        obs_new = add_noise_floor_obs_quarter(
            obs, infile=infile, scale=scale / 2.27)
    elif optype == "dime":   # constant noise floor of 10 mJy * scale
        obs_new = obs.copy()
        noise_floor = 1e-2 * scale / 2.27
        ref_noise = noise_floor * np.ones(obs_new.data['sigma'].shape)
        # add refractive noise in quadrature
        obs_new.data['sigma'] = np.sqrt(
            obs_new.data['sigma']**2 + ref_noise**2)
        obs_new.data['qsigma'] = np.sqrt(
            obs_new.data['qsigma']**2 + ref_noise**2)
        obs_new.data['usigma'] = np.sqrt(
            obs_new.data['usigma']**2 + ref_noise**2)
        obs_new.data['vsigma'] = np.sqrt(
            obs_new.data['vsigma']**2 + ref_noise**2)
    else:
        print("optype %s is not supported!" % (optype))
        obs_new = obs.copy()

    return obs_new


def add_noise_floor_obs_quarter(obs, infile, scale=1.0):
    """
    Add u,v dependent noisefloor to obs data by loading from a template
    and interpolate the noise values to match the shape of the input Obsdata.
    the templated file is in ehtim.obsdata.Obsdata format with
    the renormalized refractive noise values stored in the "vis" column

    Args:
        obs: the input obsdata
        infile: Pickle file for the model noise-floor table
        scale (float, optional):
            The scaling factor of the noise floor.
            Defaults to 1.0.

    Returns:
        obsdata: scattering-budget added obsfile.
    """
    # read obs.data and convert to a pandas dataframe
    df = pd.DataFrame(obs.data)

    # load the noise floor template
    pfile = open(infile, "rb")
    temp = pickle.load(pfile)
    nf = pd.DataFrame(temp.data)

    # interpolate the noise floor values and add to sigma for each baseline
    blgroup = df.groupby(['t1', 't2'])    # group by baseline

    df_new = pd.DataFrame()
    for key in list(blgroup.groups.keys()):
        # single baseline data
        dt = blgroup.get_group(key)
        try:
            # get the noise floor and time stamp from the template
            nfl_org = nf.groupby(['t1', 't2']).get_group(key)["vis"]
            t_org = nf.groupby(['t1', 't2']).get_group(key)["time"]
        except KeyError:
            key2 = (key[1], key[0])
            # get the noise floor and time stamp from the template
            nfl_org = nf.groupby(['t1', 't2']).get_group(key2)["vis"]
            t_org = nf.groupby(['t1', 't2']).get_group(key2)["time"]

        # interpolate
        f = interp.interp1d(t_org, nfl_org, fill_value="extrapolate")
        t_new = dt['time']
        nfl_new = f(t_new) * scale

        # add noise floor
        dt.loc[:, "sigma"] = np.sqrt(
            dt["sigma"]**2 + np.abs(nfl_new)**2) * np.sign(dt["sigma"])
        dt.loc[:, "qsigma"] = np.sqrt(
            dt["qsigma"]**2 + np.abs(nfl_new)**2) * np.sign(dt["qsigma"])
        dt.loc[:, "usigma"] = np.sqrt(
            dt["usigma"]**2 + np.abs(nfl_new)**2) * np.sign(dt["usigma"])
        dt.loc[:, "vsigma"] = np.sqrt(
            dt["vsigma"]**2 + np.abs(nfl_new)**2) * np.sign(dt["vsigma"])

        # append the new data of the baseline to the new dataframe
        df_new = df_new.append(dt)

    df_new = df_new.sort_index()    # restore the original order of the data

    # write the new sigma to data
    obs_new = obs.copy()
    obs_new.data["sigma"] = df_new["sigma"]
    obs_new.data["qsigma"] = df_new["qsigma"]
    obs_new.data["usigma"] = df_new["usigma"]
    obs_new.data["vsigma"] = df_new["vsigma"]

    return obs_new


# Reblur observations
def Blur_obs(sm, obs, use_approximate_form=True):
    """
    Blurs the observation obs by multiplying visibilities by the ensemble-average scattering kernel.
    See Fish et al. (2014): arXiv: 1409.4690.
        Args:
            sm(so.ScatteringModel): ScatteringModel
            obs(Obsdata): The observervation data(including scattering).

        Returns:
            obsblur(Obsdata): The blurred observation.
    """

    # make a copy of observation data
    datatable = (obs.copy()).data

    vis = datatable['vis']
    qvis = datatable['qvis']
    uvis = datatable['uvis']
    vvis = datatable['vvis']
    sigma = datatable['sigma']
    qsigma = datatable['qsigma']
    usigma = datatable['usigma']
    vsigma = datatable['vsigma']
    u = datatable['u']
    v = datatable['v']

    # divide visibilities by the scattering kernel
    for i in range(len(vis)):
        ker = sm.Ensemble_Average_Kernel_Visibility(
            u[i], v[i],
            wavelength_cm=299792458/obs.rf*100.0,
            use_approximate_form=use_approximate_form)
        vis[i] = vis[i] * ker
        qvis[i] = qvis[i] * ker
        uvis[i] = uvis[i] * ker
        vvis[i] = vvis[i] * ker
        sigma[i] = sigma[i] * ker
        qsigma[i] = qsigma[i] * ker
        usigma[i] = usigma[i] * ker
        vsigma[i] = vsigma[i] * ker

    datatable['vis'] = vis
    datatable['qvis'] = qvis
    datatable['uvis'] = uvis
    datatable['vvis'] = vvis
    datatable['sigma'] = sigma
    datatable['qsigma'] = qsigma
    datatable['usigma'] = usigma
    datatable['vsigma'] = vsigma

    obsblur = eh.obsdata.Obsdata(
        obs.ra, obs.dec, obs.rf, obs.bw, datatable, obs.tarr,
        source=obs.source, mjd=obs.mjd,
        ampcal=obs.ampcal, phasecal=obs.phasecal,
        opacitycal=obs.opacitycal, dcal=obs.dcal, frcal=obs.frcal)

    return obsblur
