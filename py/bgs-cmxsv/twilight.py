#!/usr/bin/env python                                                                                                                                                                                                                        
from    __future__               import  absolute_import, division, print_function

import time
import pickle
import specsim
import numpy as np

from   pkg_resources import resource_filename
from   scipy.interpolate import interp1d

import  os
import  time
import  argparse
import  fitsio
import  desisurvey
import  warnings
import  specsim.config
import  ephem
import  speclite
import  numpy                    as      np
import  astropy.units            as      u
import  pylab                    as      pl
import  matplotlib.pyplot        as      plt
import  specsim.simulator        as      simulator
import  speclite.filters         as      filters
import  astropy.units            as      units

from    astropy.table            import  Table
from    specsim                  import  config
from    scipy                    import  ndimage
from    multiprocessing          import  Pool, Array
from    desisurvey.utils         import  get_location
from    astropy.time             import  Time
from    astropy.coordinates      import  SkyCoord, EarthLocation, AltAz
from    speclite                 import  filters
from    desispec.interpolation   import  resample_flux
from    specsim.atmosphere       import  krisciunas_schaefer, Moon
from    get_sky                  import  get_sky
from    pkg_resources            import  resource_filename
from    desiutil.iers            import  freeze_iers


freeze_iers()

def XX(zd):
    # Airmass at given zenith distance.                                                                                                                                                                                                      
    return  (1. - 0.96 * np.sin(zd * np.pi / 180.)**2.)**-0.5

def get_twi(wave, alpha, delta, airmass, check=False):
    # wave: angstroms.
    # alpha: sun altitude [degrees].
    # delta: sun separation [degrees].
    # airmass: X

    assert  (alpha < -14.) & (alpha > -20.)
        
    # https://desi.lbl.gov/trac/attachment/wiki/BrightGalaxyWG/LunarModel/fagrelius_thesis_final_cochair.pdf
    ftwi       = resource_filename('bgs-cmxsv', 'dat/twilight_coeffs.p')
    twi_coeffs = pickle.load(open(ftwi, 'rb'))

    twi        = (
                  twi_coeffs['t0'] * np.abs(alpha) +      # CT2
                  twi_coeffs['t1'] * np.abs(alpha)**2 +   # CT1
                  twi_coeffs['t2'] * np.abs(delta)**2 +   # CT3
                  twi_coeffs['t3'] * np.abs(delta)        # CT4
                 ) * np.exp(-twi_coeffs['t4'] * airmass) + twi_coeffs['c0']

    twi_wave     = twi_coeffs['wave']
    twi          = np.array(twi)

    if check:
        from   feasibgs.skymodel import _cI_twi
        
        chh      = _cI_twi(-15., 90., 1.3)[1]

        assert  np.allclose(chh, twi)
    
    twi         /= np.pi

    # fill_value='extrapolate'.
    I_twi_interp = interp1d(10. * twi_wave, twi, bounds_error=False, fill_value='extrapolate')

    # Clip below. 
    Isky         = np.clip(I_twi_interp(wave), 0.0, None)

    # [$10^{-17} erg/cm^{2}/s/\AA/arcsec^2$]
    return  twi_wave, Isky


if __name__ == '__main__':
    import pylab as pl

    start             = time.time()

    # specsim.
    config            = specsim.config.load_config('desi')
    simulator         = specsim.simulator.Simulator('desi')
    
    simulator.simulate()

    rfilter           = filters.load_filters('decam2014-r')
    
    gfas              = resource_filename('bgs-cmxsv', 'dat/offline_all_guide_ccds_SV1-thru_20210111.fits')
    gfas              = Table.read(gfas)

    # 
    # gfas            = gfas[gfas['GFA_TRANSPARENCY_MED'] > 0.95] 

    mayall            = get_location()

    emayall           = ephem.Observer()
    emayall.lon       = ephem.degrees(mayall.lon.value * np.pi / 180.)
    emayall.lat       = ephem.degrees(mayall.lat.value * np.pi / 180.)
    emayall.elevation = mayall.height.value

    moon              = ephem.Moon()
    sun               = ephem.Sun()

    fullwave          = config.wavelength

    '''
    fig, axes         = plt.subplots(1, 3, figsize=(15, 5))
    
    xs                = np.arange(17., 21., 0.1)

    axes[0].set_ylabel('Sun sep. [deg.]')
    
    for ax in axes:
        ax.plot(xs, xs, c='k', lw=0.2)

        ax.set_xlim(-14., -20.0)
        ax.set_ylim(0.0,   180.)

        ax.set_xlabel('Sun alt. [deg.]')
    ''' 
    print('\n\n')

    print('MJD \t\t NIGHT \t\t EXPID \t\t CAMERA \t EXPTIME \t SUNALT \t SUNSEP \t ZD \t\t X \t\t TRANS \t\t GFA_R \t\t MODEL_R \t GFA_R - MODEL_R')

    first_night = None
    
    for	i, (night, mjd, expid, camera, exptime, ra, dec, trans, gfa_r) in enumerate(zip(gfas['NIGHT'], gfas['MJD'], gfas['EXPID'], gfas['EXTNAME'], gfas['EXPTIME'],\
                                                                                        gfas['RACEN'], gfas['DECCEN'], gfas['TRANSPARENCY'], gfas['SKY_MAG_AB'])):        

        if first_night == None:
            first_night = night

        t                 = Time(mjd, format='mjd', scale='utc')
        emayall.date      = t.iso

        pos               = SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame='icrs').transform_to(AltAz(obstime=t, location=mayall))
        alt               = pos.alt.degree
        az                = pos.az.degree
        zd                = 90. - alt

        sun.compute(emayall)
        moon.compute(emayall)

        sun_alt           = sun.alt * (180. / np.pi)
        sun_sep           = desisurvey.utils.separation_matrix([ra] * u.deg, [dec] * u.deg, [sun.ra] * u.deg, [sun.dec] * u.deg)[0][0].value

        if (sun_alt > -14.) | (sun_alt < -20.):
            continue
        
        moon_alt          = moon.alt * (180. / np.pi)
        moon_frac         = moon.moon_phase
        moon_sep          = desisurvey.utils.separation_matrix([ra] * u.deg, [dec] * u.deg, [moon.ra] * u.deg, [moon.dec] * u.deg)[0][0].value
        moon_phase        = 180. * np.arccos(2.*moon_frac - 1.) / np.pi # deg.                                                                                                                                                                    
        moon_zd           = (90. - moon_alt)

        X                 = XX(zd)
        '''
        # 
        simulator.atmosphere.airmass                         = X
        simulator.atmosphere.moon.moon_zenith                = moon_zd * u.deg
        simulator.atmosphere.moon.separation_angle           = moon_sep * u.deg
        simulator.atmosphere.moon.moon_phase                 = moon_phase / 180.

        simulator.simulate()

        notwi_sky         = simulator.atmosphere.surface_brightness
        notwi_sky        *= 1.e-17

        twiwave, twi      = get_twi(fullwave, sun_alt, sun_sep, X, check=False)
        twi              *= 1.e-17

        sky               = twi + notwi_sky.value
        sky              *= u.erg / (u.cm ** 2 * u.s * u.angstrom)
        
        sky_pad, skywave_pad   = rfilter.pad_spectrum(sky, fullwave.value, method="zero")

        # Normalize to Dark Sky Zenith V [u.erg / (u.cm ** 2 * u.s * u.angstrom]                                                                                                                                                            
        model_rmag             = rfilter.get_ab_magnitudes(sky_pad, skywave_pad).as_array()[0][0]
        '''

        model_rmag = 99.
        
        print('{:.6f} \t {} \t {:08d} \t {} \t {:.6f} \t {:.6f} \t {:.6f} \t {:.6f} \t {:.6f} \t {:.6f} \t {:.6f} \t {:.6f} \t {:.6f}'.format(mjd, night, expid, camera, exptime, sun_alt, sun_sep,\
                                                                                                                                              zd, X, trans, gfa_r, model_rmag, gfa_r - model_rmag))

        # zero = axes[0].scatter(sun_alt, sun_sep, c=gfa_r,              marker='.', s=12, vmin=18., vmax=21.)
        # one  = axes[1].scatter(sun_alt, sun_sep, c=gfa_r - model_rmag, marker='.', s=12, vmin=-.5, vmax=.5)
        # two  = axes[2].scatter(sun_alt, sun_sep, c=X,                  marker='.', s=12, vmin=1.,  vmax=2.)

        if night != first_night:
            exit(0)
'''        
fig.colorbar(zero, ax=axes[0], label='GFA R')
fig.colorbar( one, ax=axes[1], label='GFA R - MODEL R')
fig.colorbar( two, ax=axes[2], label='AIRMASS')

pl.show()
'''
print('\n\n')
