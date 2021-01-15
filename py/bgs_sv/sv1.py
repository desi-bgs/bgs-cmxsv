'''


module for read in useful SV1 data 


'''
import os 
import yaml
import fitsio
import numpy as np 
# -- astropy -- 
from astropy.io import fits 
import astropy.table as atable
import astropy.units as units
import astropy.constants as constants


assert os.environ['NERSC_HOST'] == 'cori'
assert os.environ['DESIMODEL'] != ''


def sv1_exposures(): 
    ''' read in  summary of SV1 exposures (Anand, from Aaron & David K) to
    astropy table. `sv1-exposures.fits` was compiled by Mike W. 
    '''
    fexp = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat', 'sv1-exposures.fits')
    return atable.Table.read(fexp)


def sv1_deep_exposures(): 
    ''' read in summary of SV1 *deep* exposures to astropy table. These
    deep "exposures" are combination of a bunch of separate exposures. 
    '''
    # read all sv1 exposures 
    sv1 = sv1_exposures() 

    # read in tileid, date, exposure id
    _, _, deep_expid = np.loadtxt(
            os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat',
                'blanc_deep_explist.dat'), 
            unpack=True)
    # match them to the sv1 exposure table using the expid 
    deep = sv1[np.isin(sv1['EXPID'], deep_expid)] 
   
    # tile information 
    tinfo = atable.unique(deep['TILEID', 'TILERA', 'TILEDEC'], keys='TILEID')

    # bin exposures by tileid and sum up the exposure time and depths
    dexps = deep['TILEID', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']
    dexps_binned = dexps.group_by('TILEID').groups.aggregate(np.sum)
    
    return atable.join(dexps_binned, tinfo, keys='TILEID', join_type='left')


def sv1_bright_exposures(): 
    ''' read in only the bright time SV1 exposures.
    '''
    exps = sv1_exposures()
    
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
    dir_thru = os.path.join(os.getenv("DESIMODEL"), 'data', 'throughput') 
     
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
        norm_cam    = np.zeros(len(waves[camera]))
        sky_cam     = np.zeros(len(waves[camera]))

        for spec in specs:
            if ftype == "data":
                raise NotImplementedError
            elif ftype == "model":
                fn = os.path.join(reduxdir, "exposures", night, expid, 
                        "sky-%s%s-%s.fits" % (camera, spec, expid),)
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


def get_zbest(tileid, date, expid=None, redux='blanc', targetclass='all'): 
    ''' read good fibers that point to BGS targets in zbest files of (tileid,
    date). 

    Notes
    -----
    * for deep exposure set date = 'deep'
    * at the moment it only support the nightly coadds so expid=None 
    '''
    from desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
    from desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask

    dir_redux = '/global/cfs/cdirs/desi/spectro/redux/%s/' % redux
    dir_zbest = os.path.join(dir_redux, 'tiles', str(tileid), str(date)) 
    
    petals = [] 
    for petal in range(10):
        fzbest = os.path.join(dir_zbest, 'zbest-%i-%i-%s.fits' % (petal,
            tileid, str(date))) 

        if not os.path.isfile(fzbest): 
            print(" %s does not exist" % fzbest) 
            continue 

        infile  = fits.open(fzbest)
        zbest   = atable.Table(infile['ZBEST'].data) 
        fmap    = atable.Table(infile['FIBERMAP'].data) 
        
        tinfo = fmap['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'FLUX_R',
                'FIBERFLUX_R', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET',
                'DESI_TARGET', 'BGS_TARGET'] 
        tinfo    = atable.unique(tinfo, keys='TARGETID')
        
        deep     = atable.join(zbest, tinfo, keys='TARGETID', join_type='left')
        assert len(deep) == 500 # fibers per petal
        
        # bad fibers
        badfiber = (deep['ZWARN'] & 2**9) != 0 # no data 

        # only keep bgs 
        badbgs   = (deep['SV1_DESI_TARGET'] & sv1_desi_mask['BGS_ANY']) == 0

        nbright  = (deep['SV1_BGS_TARGET'] & sv1_bgs_mask['BGS_BRIGHT']) == 0
        nfaint   = (deep['SV1_BGS_TARGET'] & sv1_bgs_mask['BGS_FAINT'])  == 0
        
        if targetclass == 'all':
            # limit to faint | bright bgs only.
            badbgs   = badbgs | (nbright & nfaint)
        elif targetclass == 'bright': 
            # limit to bright bgs only.
            badbgs   = badbgs | nbright
        else: 
            raise NotImplementedError
    
        cuts = ~badfiber & ~badbgs
        petals.append(deep[~cuts])
    
    return atable.vstack(petals) 
