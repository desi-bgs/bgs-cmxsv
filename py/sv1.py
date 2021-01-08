'''


module for read in useful SV1 data 


'''
import os 
import yaml
import numpy as np 
# -- astropy -- 
from astropy.table import Table
import astropy.constants as constants


assert os.environ['NERSC_HOST'] == 'cori'


def sv1_exposures(): 
    ''' read in  summary of SV1 exposures (Anand, from Aaron & David K) to
    astropy table. `/global/homes/m/mjwilson/desi/SV1/sv1-exposures.fits` was
    compiled by Mike W. 
    '''
    fexp = '/global/homes/m/mjwilson/desi/SV1/sv1-exposures.fits'

    return Table.read('sv1-exposures.fits')


def sv1_bright_exposures(): 
    ''' read in only the bright time SV1 exposures.
    '''
    exps = sv1_exposure()
    
    # calculate bright limit 
    nominal_dark = 21.07
    dmab = -2.5 * np.log10(2.5)
    bright_lim = nominal_dark + dmab
    
    is_bright = (exps['GFA_SKY_MAG_AB_MED'] <= bright_lim) 
    return exps[is_bright]


def get_obs_sky(night, expid, exptime, ftype="model", redux="daily", smoothing=100.0):
    ''' get observed sky surface brightness for exposure specified by (night,
    expid, exptime) 
    
    # AR r-band sky mag / arcsec2 from sky-....fits files
    # AR ftype = "data" or "model"
    # AR redux = "daily" or "blanc"
    # AR if ftype = "data" : read the sky fibers from frame*fits + apply flat-field
    # AR if ftype = "model": read the sky model from sky*.fits for the first fiber of each petal (see DJS email from 29Dec2020)
    # AR those are in electron / angstrom
    # AR to integrate over the decam-r-band, we need cameras b and r
    '''
    assert ftype in ['data', 'model'], "ftype should be 'data' or 'model'" 
    night = str(night)
    expid = str(expid).zfill(8)


    # AR/DK DESI spectra wavelengths
    wmin, wmax, wdelta = 3600, 9824, 0.8
    wave_full = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
    cslice  = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}
    waves   = {'b': wave_full[cslice['b']], 'r': wave_full[cslice['r']], 'z': wave_full[cslice['z']]}

    # AR DESI telescope geometric area (cm2) and fiber area (arcsec2)
    f   = open(os.path.join(os.getenv("DESIMODEL"), "data", "desi.yaml"), "r")
    desi = yaml.safe_load(f)
    f.close()
    telap_cm2 = desi["area"]["geometric_area"] * 1e4  # AR telescope geometric area in cm2
    fiber_area_arcsec2 = (np.pi * (desi["fibers"]["diameter_arcsec"] / 2.0) ** 2)  # fiber

    # load spectrograph throughput and interpolate to wavelength 
    dir_thru = os.path.join(os.environ['DESIMODEL'], '/data/throughput/') 
     
    thru = {}
    for camera in ['b', 'r', 'z']: 
        data = fitsio.read(os.path.join(dir_thru, 'thru-%s.fits' % camera), 'THROUGHPUT') 
        thru[camera] = np.interp(waves[camera], data["wavelength"], data["throughput"])

    
    dailydir = "/global/cfs/cdirs/desi/spectro/redux/daily/"
    reduxdir    = dailydir.replace("daily", redux)
    
    # AR see AK email [desi-data 5218]
    sky = np.zeros(len(wave_full))
    if redux == "blanc":
        specs = ["0", "1", "3", "4", "5", "7", "8", "9"]
    else:
        specs = np.arange(10, dtype=int).astype(str)

    for camera in ["b", "r", "z"]:
        norm_cam    = np.zeros(len(waves[camera]]))
        sky_cam     = np.zeros(len(waves[camera]]))

        for spec in specs:
            if ftype == "data":
                raise NotImplementedError
            elif ftype == "model":
                fn = os.path.join(reduxdir, "exposures", night, expid, 
                        "sky-%s%i-%s.fits" % (camera, spec, expid),)
                if not os.path.isfile(fn):
                    print("Skipping non-existent %s" % fn)
                    continue 

                fd = fitsio.FITS(fn)
                with fd as hdus:
                    exptime = hdus[0].read_header()['EXPTIME']

                    flux = hdus['SKY'].read()
                    ivar = hdus['IVAR'].read()
                    mask = hdus['MASK'].read()
                    
                    # Verify that we have the expected wavelengths.
                    assert np.allclose(waves[camera], fd["WAVELENGTH"].read())

                    # Verify that ivar is purely statistical.
                    # assert np.array_equal(ivar, hdus['STATIVAR'].read())
                    # Verify that the model has no masked pixels.
                    # assert np.all((mask == 0) & (ivar > 0))
                    # Verify that the sky model is constant.
                    # assert np.array_equal(np.max(ivar, axis=0), np.min(ivar, axis=0))

                    # assert np.allclose(np.max(flux, axis=0), np.min(flux, axis=0))
                    # There are actually small variations in flux!
                    # TODO: figure out where these variations come from.
                    #if fd["IVAR"][0, :][0].max() > 0:
                    #    sky_cam  += fd["SKY"][0, :][0]  # AR reading the first fiber only
                    #    norm_cam += np.ones(len(fullwave[cslice[camera]]))
                    
                    # For now, take the median over fibers.
                    sky_cam  += np.median(flux, axis=0)
                    norm_cam += np.ones(len(waves[camera]))

        # AR sky model flux in incident photon / angstrom / s
        # if nspec > 0:
        keep = norm_cam > 0

        if keep.sum() == 0:
            print("%s-%s-%s: no spectra for %s" % (night, expid, camera, ftype))
        sky[cslice[camera]][keep] = (sky_cam[keep] / norm_cam[keep] / exptime / thru[camera][keep])
    
    # convert sky to erg/A/s/cm^2/arcsec^2

    # AR sky model flux in erg / angstrom / s (using the photon energy in erg)
    e_phot_erg = (
        constants.h.to(units.erg * units.s)
        * constants.c.to(units.angstrom / units.s)
        / (wave_full * units.angstrom)
    )
    sky *= e_phot_erg.value
    # AR sky model flux in erg / angstrom / s / cm**2 / arcsec**2
    sky /= telap_cm2 * fiber_area_arcsec2
    return wave_full, sky
