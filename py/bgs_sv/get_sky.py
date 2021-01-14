import os
import yaml
import glob 
import ephem
import fitsio
import pickle
import desisurvey
import astropy.constants as constants
import astropy.units as units
import numpy as np 
import pylab as pl
import scipy

from desisurvey.utils import get_location, get_date
from astropy.coordinates import Angle
from astropy.time import Time
from desispec.fiberflat import apply_fiberflat

from scipy.interpolate import interp1d
from speclite import filters
from multiprocessing import Pool
from functools import partial
from desispec.io.fluxcalibration import read_average_flux_calibration
from desispec.calibfinder import  CalibFinder
from desispec.io import  read_frame, read_fiberflat, read_flux_calibration, read_sky, read_fibermap

# -- astropy --
import astropy.time
import astropy.coordinates

from astropy.io import fits 
from astropy import units as u
from astropy.table import Table

# -- eso sky --
# from skycalc_cli import skycalc 

# -- feasibgs -- 
from feasibgs import util as UT
from feasibgs import skymodel as Sky
from feasibgs.skymodel import Isky_newKS_twi, _specsim_initialize

# -- plotting -- 
import matplotlib as mpl
import matplotlib.pyplot as plt

# -- specsim --
import specsim
import specsim.simulator as      simulator

from astropy.table import  Table
from specsim import  config

# -- desispec --
from desispec.interpolation import resample_flux


config = specsim.config.load_config('desi')
simulator = specsim.simulator.Simulator('desi')

simulator.simulate()

# AR DESI telescope geometric area (cm2) and fiber area (arcsec2)
# AR for computing SKY_RMAG_AB
fn = os.path.join(os.getenv("DESIMODEL"), "data", "desi.yaml")
f = open(fn, "r")
desi = yaml.safe_load(f)
f.close()
telap_cm2 = desi["area"]["geometric_area"] * 1e4  # AR telescope geometric area in cm2
fiber_area_arcsec2 = (
    np.pi * (desi["fibers"]["diameter_arcsec"] / 2.0) ** 2
)  # fiber area in arcsec2

# AR/DK exposure depths utilities
def load_spec_thru(path=os.getenv("DESIMODEL") + "/data/throughput/"):
    thru = {}
    for camera in ["b", "r", "z"]:
        data = fitsio.read(f"{path}/thru-{camera}.fits", "THROUGHPUT")
        thru[camera] = np.interp(
            fullwave[cslice[camera]], data["wavelength"], data["throughput"]
        )
    return thru

# AR folders / files
sv1dir = "/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/"
dailydir = "/global/cfs/cdirs/desi/spectro/redux/daily/"
pixwfn = "/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/pixweight/sv1/resolve/dark/sv1pixweight-dark.fits"
desfn = os.path.join(sv1dir, "misc", "des_footprint.txt")
gfafn = np.sort(
    glob.glob(
        "/global/cfs/cdirs/desi/users/ameisner/GFA/conditions/offline_all_guide_ccds_SV1-thru_20??????.fits"
    )
)[-1]

# AR/DK DESI spectra wavelengths
wmin, wmax, wdelta = 3600, 9824, 0.8
fullwave = np.round(np.arange(wmin, wmax + wdelta, wdelta), 1)
cslice = {"b": slice(0, 2751), "r": slice(2700, 5026), "z": slice(4900, 7781)}

spec_thru = load_spec_thru()

nightly_dsky_cache = {}

# AR r-band sky mag / arcsec2 from sky-....fits files
def get_sky(night, expid, exptime, ftype="model", redux="daily", smoothing=100.0, specsim_darksky=False, nightly_darksky=False):
    # AR ftype = "data" or "model"
    # AR redux = "daily" or "blanc"
    # AR if ftype = "data" : read the sky fibers from frame*fits + apply flat-field
    # AR if ftype = "model": read the sky model from sky*.fits for the first fiber of each petal (see DJS email from 29Dec2020)
    # AR those are in electron / angstrom
    # AR to integrate over the decam-r-band, we need cameras b and r
    if ftype not in ["data", "model"]:
        sys.exit("ftype should be 'data' or 'model'")

    sky = np.zeros(len(fullwave))
    reduxdir = dailydir.replace("daily", redux)
    
    # AR see AK email [desi-data 5218]
    if redux == "blanc":
        specs = ["0", "1", "3", "4", "5", "7", "8", "9"]
    else:
        specs = np.arange(10, dtype=int).astype(str)
        
    for camera in ["b", "r", "z"]:
        norm_cam = np.zeros(len(fullwave[cslice[camera]]))
        sky_cam = np.zeros(len(fullwave[cslice[camera]]))
        
        for spec in specs:
            # AR data
            if ftype == "data":
                frfn = os.path.join(
                    reduxdir,
                    "exposures",
                    "{}".format(night),
                    expid,
                    "frame-{}{}-{}.fits".format(camera, spec, expid),
                )

                if not os.path.isfile(frfn):
                    print("Skipping non-existent {}".format(frfn))

                    continue
                
                fr = read_frame(frfn, skip_resolution=True)
                    
                flfn = os.environ['DESI_SPECTRO_CALIB'] + '/' + CalibFinder([fr.meta]).data['FIBERFLAT']
 
                # calib_flux [1e-17 erg/s/cm2/Angstrom] = uncalib_flux [electron/Angstrom] / (calibration_model * exptime [s])
                # fl = read_fiberflat(flfn)
 
                # No exptime correction. 
                # apply_fiberflat(fr, fl)
                        
                # AR cutting on sky fibers with at least one valid pixel.
                ii = (fr.fibermap["OBJTYPE"] == "SKY") & (fr.ivar.sum(axis=1) > 0)
                
                # AR frame*fits are in e- / angstrom ; adding the N sky fibers
                # sky_cam += fr.flux[ii, :].sum(axis=0)
                # nspec += ii.sum()
                
                # Ignores fiberflat corrected (fl), as includes e.g. fiberloss. 
                sky_cam  += (fr.flux[ii, :] * fr.ivar[ii, :]).sum(axis=0)
                norm_cam += fr.ivar[ii, :].sum(axis=0)
                    
            # AR model
            if ftype == "model":
                fn = os.path.join(
                    reduxdir,
                    "exposures",
                    "{}".format(night),
                    expid,
                    "sky-{}{}-{}.fits".format(camera, spec, expid),
                )
                if not os.path.isfile(fn):
                    print("Skipping non-existent {}".format(fn))
                    
                else:
                    print("Solving for {}".format(fn))
                    
                    fd = fitsio.FITS(fn)
                    
                    assert np.allclose(
                        fullwave[cslice[camera]], fd["WAVELENGTH"].read()
                    )
                    
                    fd = fitsio.FITS(fn)

                    with fd as hdus:
                        exptime = hdus[0].read_header()['EXPTIME']
                        
                        flux = hdus['SKY'].read()
                        ivar = hdus['IVAR'].read()
                        mask = hdus['MASK'].read()
                
                        # Verify that we have the expected wavelengths.
                        # assert np.allclose(detected[camera].wave, hdus['WAVELENGTH'].read())
                        # Verify that ivar is purely statistical.
                        # assert np.array_equal(ivar, hdus['STATIVAR'].read())
                        # Verify that the model has no masked pixels.
                        # assert np.all((mask == 0) & (ivar > 0))
                        # Verify that the sky model is constant.
                        # assert np.array_equal(np.max(ivar, axis=0), np.min(ivar, axis=0))
              
                        # assert np.allclose(np.max(flux, axis=0), np.min(flux, axis=0))
                        # There are actually small variations in flux!
                        # TODO: figure out where these variations come from.
                        # For now, take the median over fibers.
                        if fd["IVAR"][0, :][0].max() > 0:
                            sky_cam  += fd["SKY"][0, :][0]  # AR reading the first fiber only
                            norm_cam += np.ones(len(fullwave[cslice[camera]]))
                        
                        # sky_cam  += np.median(flux, axis=0)
                        # norm_cam += np.ones(len(fullwave[cslice[camera]]))
                    
                    fd.close()

        # AR sky model flux in incident photon / angstrom / s
        # if nspec > 0:
        keep = norm_cam > 0
        
        if keep.sum() > 0:
            # [e/A/s] / throughput.
            sky[cslice[camera]][keep] = (
                sky_cam[keep] / norm_cam[keep] / exptime / spec_thru[camera][keep]
            )
        else:
            print("{}-{}-{}: no spectra for {}".format(night, expid, camera, ftype))
            
    # AR sky model flux in erg / angstrom / s (using the photon energy in erg).
    e_phot_erg = (
        constants.h.to(units.erg * units.s)
        * constants.c.to(units.angstrom / units.s)
        / (fullwave * units.angstrom)
    )
    sky *= e_phot_erg.value
    
    # AR sky model flux in [erg / angstrom / s / cm**2 / arcsec**2].
    sky /= (telap_cm2 * fiber_area_arcsec2)
    
    if smoothing > 0.0: 
        sky = scipy.ndimage.gaussian_filter1d(sky, 100.) 
        
    # AR integrate over the DECam r-band
    vfilter = filters.load_filters('bessell-V')
    rfilter = filters.load_filters('decam2014-r')
        
    # AR zero-padding spectrum so that it covers the DECam r-band range
    sky_pad, fullwave_pad = vfilter.pad_spectrum(sky, fullwave, method="zero")
    vmag    = vfilter.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), fullwave_pad * units.angstrom).as_array()[0][0]

    # AR zero-padding spectrum so that it covers the DECam r-band range
    sky_pad, fullwave_pad = rfilter.pad_spectrum(sky, fullwave, method="zero")
    rmag    = rfilter.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), fullwave_pad * units.angstrom).as_array()[0][0]

    if specsim_darksky:
        fd  = fitsio.FITS(fn)
        
        # Dark Sky at given airmass (cnst. across spectrograph / camera).
        simulator.atmosphere.airmass = fd['SKY'].read_header()['AIRMASS']

        # Force below the horizon: dark sky.                                                                                                                                                                                               
        simulator.atmosphere.moon.moon_zenith = 120. * u.deg

        simulator.simulate()

        # [1e-17 erg / (Angstrom arcsec2 cm2 s)].
        sim_darksky  = simulator.atmosphere.surface_brightness
        sim_darksky *= 1.e-17

        dsky_pad, dskywave_pad = rfilter.pad_spectrum(sim_darksky.value, config.wavelength.value, method="zero")
        dsky_rmag              = rfilter.get_ab_magnitudes(dsky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), dskywave_pad * units.angstrom).as_array()[0][0]
        
        sim_darksky  = resample_flux(fullwave, config.wavelength.value, sim_darksky.value)
        sim_darksky *= u.erg / (u.cm ** 2 * u.s * u.angstrom * u.arcsec ** 2)
        
        sky          = np.clip(sky - sim_darksky.value, a_min=0.0, a_max=None)

        # AR zero-padding spectrum so that it covers the DECam r-band range
        sky_pad, fullwave_pad = vfilter.pad_spectrum(sky, fullwave, method="zero")
        vmag_nodark  = vfilter.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), fullwave_pad * units.angstrom).as_array()[0][0]

        # AR zero-padding spectrum so that it covers the DECam r-band range                                                                                                                                                                
        sky_pad, fullwave_pad = rfilter.pad_spectrum(sky, fullwave, method="zero")
        rmag_nodark  = rfilter.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), fullwave_pad * units.angstrom).as_array()[0][0]

        return  fullwave, sky, vmag, rmag, vmag_nodark, rmag_nodark

    elif nightly_darksky:
        from   pkg_resources import resource_filename

        
        if night in nightly_dsky_cache.keys():
            return  nightly_dsky_cache[night]
        
        gfa_info       = resource_filename('bgs-cmxsv', 'dat/sv1-exposures.fits')
        gfa_info       = Table.read(gfa_info)

        bright_cut     = 20.50
        
        good_conds     = (gfa_info['GFA_TRANSPARENCY_MED'] > 0.95) & (gfa_info['GFA_SKY_MAG_AB_MED'] >= bright_cut)
        good_conds     =  gfa_info[good_conds]
        
        expids         = np.array([np.int(x.split('/')[-1]) for x in glob.glob(os.path.join(reduxdir, "exposures", "{}/00*".format(night)))])
        isgood         = np.isin(good_conds['EXPID'], expids)

        if np.any(isgood):        
            good_conds = good_conds[isgood]
            best       = good_conds['GFA_SKY_MAG_AB_MED'] == good_conds['GFA_SKY_MAG_AB_MED'].max()

            print('Nightly Dark GFA r-mag for {}:  {}'.format(night, good_conds['GFA_SKY_MAG_AB_MED'].max()))
            
            best_expid = good_conds[best]['EXPID'][0]
            best_expid = '{:08d}'.format(best_expid)
            
            darkwave, darksky, dark_vmag, dark_rmag = get_sky(night, best_expid, exptime, ftype="model", redux="daily", smoothing=0.0, specsim_darksky=False, nightly_darksky=False)

            sky        = np.clip(sky - darksky, a_min=0.0, a_max=None)

            # AR zero-padding spectrum so that it covers the DECam r-band range
            sky_pad, fullwave_pad = vfilter.pad_spectrum(sky, fullwave, method="zero")
            vmag_nodark  = vfilter.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), fullwave_pad * units.angstrom).as_array()[0][0]

            # AR zero-padding spectrum so that it covers the DECam r-band range
            sky_pad, fullwave_pad = rfilter.pad_spectrum(sky, fullwave, method="zero")
            rmag_nodark  = rfilter.get_ab_magnitudes(sky_pad * units.erg / (units.cm ** 2 * units.s * units.angstrom), fullwave_pad * units.angstrom).as_array()[0][0]
            
            nightly_dsky_cache[night] = fullwave, sky, vmag, rmag, vmag_nodark, rmag_nodark
            
        else:
            print('No nightly darksky available for night: {}.  Defaulting to specsim.'.format(night))

            nightly_dsky_cache[night] = get_sky(night, expid, exptime, ftype="model", redux="daily", smoothing=0.0, specsim_darksky=True, nightly_darksky=False)

        return  nightly_dsky_cache[night]
            
    else:
        pass
    
    return  fullwave, sky, vmag, rmag


if __name__ == '__main__':
    get_sky('20210103', '00070758', 300.094, 'model', redux="daily", smoothing=0., specsim_darksky=False, nightly_darksky=True)
    
