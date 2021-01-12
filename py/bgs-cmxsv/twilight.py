import pickle
import numpy as np

from   pkg_resources import resource_filename
from   scipy.interpolate import interp1d


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
    I_twi_interp = interp1d(10. * twi_wave, twi, bounds_error=True)

    # Clip below. 
    Isky         = np.clip(I_twi_interp(wave), 0.0, None)

    # [$10^{-17} erg/cm^{2}/s/\AA/arcsec^2$]
    return  twi_wave, Isky


if __name__ == '__main__':
    import pylab as pl
    

    wave    = np.arange(4.e3, 9.e3, 10.)
    _, twi  = get_twi(wave, -15., 90., 1.3)
        
    pl.plot(wave, twi)
    pl.show()
