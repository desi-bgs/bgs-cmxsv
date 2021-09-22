'''

module for interfacing with SV data 

'''
import os 
import glob
import yaml
import json
import fitsio
import numpy as np 
from itertools import chain, combinations_with_replacement
# -- astropy -- 
from astropy.io import fits 
import astropy.table as atable
import astropy.units as units
import astropy.constants as constants


def info_exposures(survey=None, release='everest'): 
    ''' read in exposure info in specified release 
    '''
    if release == 'everest': 
        fexp = '/global/cfs/cdirs/desi/spectro/redux/everest/exposures-everest.fits'
    else: 
        raise ValueError

    exps = atable.Table.read(fexp, hdu=1)

    if survey is None: 
        return exps
    else: 
        # only return exposures from specific surveys
        return exps[exps['SURVEY'] == survey]


def info_tiles(survey=None, release='everest'): 
    ''' read in tile info in specified release 
    '''
    if release == 'everest': 
        ftiles = '/global/cfs/cdirs/desi/spectro/redux/everest/tiles-everest.csv'
    else: 
        raise ValueError

    tiles = atable.Table.read(ftiles)

    if survey is None: 
        return tiles 
    else: 
        # only return exposures from specific surveys
        return tiles[tiles['SURVEY'] == survey]
    
def extinction_total_to_selective_ratio(band, photsys, match_legacy_surveys=False) :
    """Return the linear coefficient R_X = A(X)/E(B-V) where
    A(X) = -2.5*log10(transmission in X band),
    for band X in 'G','R' or 'Z' when
    photsys = 'N' or 'S' specifies the survey (BASS+MZLS or DECALS),
    or for band X in 'G', 'BP', 'RP' when photsys = 'G' (when gaia dr2)
    or for band X in 'W1', 'W2', 'W3', 'W4' when photsys is either 'N' or 'S'
    E(B-V) is interpreted as SFD.
    Args:
        band : 'G', 'R', 'Z', 'BP', 'RP', 'W1', 'W2', 'W3', or 'W4'
        photsys : 'N' or 'S'
    Returns:
        scalar, total extinction A(band) = -2.5*log10(transmission(band))
    """
    if match_legacy_surveys :
        # Based on the fit from the columns MW_TRANSMISSION_X and EBV
        # for the DR8 target catalogs and propagated in fibermaps
        # R_X = -2.5*log10(MW_TRANSMISSION_X) / EBV
        # It is the same value for the N and S surveys in DR8 and DR9 catalogs.
        R={"G_N":3.2140,
           "R_N":2.1650,
           "Z_N":1.2110,
           "G_S":3.2140,
           "R_S":2.1650,
           "Z_S":1.2110,
           "G_G":2.512,
           "BP_G":3.143,
           "RP_G":1.663,
        }
    else :
        # From https://desi.lbl.gov/trac/wiki/ImagingStandardBandpass
        # DECam u  3881.6   3.994
        # DECam g  4830.8   3.212
        # DECam r  6409.0   2.164
        # DECam i  7787.5   1.591
        # DECam z  9142.7   1.211
        # DECam Y  9854.5   1.063
        # BASS g  4772.1   3.258
        # BASS r  6383.6   2.176
        # MzLS z  9185.1   1.199
        # Consistent with the synthetic magnitudes and function dust_transmission

        R={"G_N":3.258,
           "R_N":2.176,
           "Z_N":1.199,
           "G_S":3.212,
           "R_S":2.164,
           "Z_S":1.211,
           "G_G":2.197,
           "BP_G":2.844,
           "RP_G":1.622,
        }

    # Add WISE from
    # https://github.com/dstndstn/tractor/blob/main/tractor/sfd.py#L23-L35
    R.update({
        'W1_N': 0.184,
        'W2_N': 0.113,
        'W3_N': 0.0241,
        'W4_N': 0.00910,
        'W1_S': 0.184,
        'W2_S': 0.113,
        'W3_S': 0.0241,
        'W4_S': 0.00910
        })

    assert(band.upper() in ["G","R","Z","BP","RP",'W1','W2','W3','W4'])
    assert(photsys.upper() in ["N","S","G"])
    return R["{}_{}".format(band.upper(),photsys.upper())]

def mwdust_transmission(ebv, band, photsys, match_legacy_surveys=False):
    """Convert SFD E(B-V) value to dust transmission 0-1 for band and photsys
    Args:
        ebv (float or array-like): SFD E(B-V) value(s)
        band (str): 'G', 'R', 'Z', 'W1', 'W2', 'W3', or 'W4'
        photsys (str or array of str): 'N' or 'S' imaging surveys photo system
    Returns:
        scalar or array (same as ebv input), Milky Way dust transmission 0-1
    If `photsys` is an array, `ebv` must also be array of same length.
    However, `ebv` can be an array with a str `photsys`.
    Also see `dust_transmission` which returns transmission vs input wavelength
    """
    if isinstance(photsys, str):
        r_band = extinction_total_to_selective_ratio(band, photsys, match_legacy_surveys=match_legacy_surveys)
        a_band = r_band * ebv
        transmission = 10**(-a_band / 2.5)
        return transmission
    else:
        photsys = np.asarray(photsys)
        if np.isscalar(ebv):
            raise ValueError('array photsys requires array ebv')
        if len(ebv) != len(photsys):
            raise ValueError('len(ebv) {} != len(photsys) {}'.format(
                len(ebv), len(photsys)))

        transmission = np.zeros(len(ebv))
        for p in np.unique(photsys):
            ii = (photsys == p)
            r_band = extinction_total_to_selective_ratio(band, p, match_legacy_surveys=match_legacy_surveys)
            a_band = r_band * ebv[ii]
            transmission[ii] = 10**(-a_band / 2.5)

        return transmission

def rr_exposure(tileid, expid, release='everest'):
    ''' redrock redshift success for given tile ID and exposure ID
    '''
    if release == 'everest': 
        dat_dir = '/global/cfs/cdirs/desi/spectro/redux/everest/tiles/perexp/'
    else: 
        raise ValueError
    
    exp_dir = os.path.join(dat_dir, str(tileid), str(expid).zfill(8))
    assert os.path.isdir(exp_dir), "%s doesn't exist" % exp_dir
    fRR = lambda i: os.path.join(exp_dir, 
            'redrock-%i-%i-exp%s.fits' % (i, tileid, str(expid).zfill(8)))
    
    petals = [] 
    for i in range(10): 
        frr = fRR(i)

        if not os.path.isfile(frr): 
            print("  %s does not exist" % frr) 
            continue 

        hdu     = fits.open(frr)
        zbest   = atable.Table(hdu['REDSHIFTS'].data) # redrock output 
        fmap    = atable.Table(hdu['FIBERMAP'].data)  # fibermap data 
        tsnr2   = atable.Table(hdu['TSNR2'].data)     # TSNR2 data 
        
        # zbest and fmap should be row-matched starting from Everest :) 
        petal   = atable.join(zbest, fmap, keys='TARGETID', join_type='left')
        petal   = atable.join(petal, tsnr2, keys='TARGETID', join_type='left')
        
        # bad fibers
        badfiber  = (petal['ZWARN'] & 2**9)  != 0 # no data 
        badfiber |= (petal['ZWARN'] & 2**11) != 0 # poor positioning
        badfiber |= (petal['ZWARN'] & 2**7)  != 0 # broken
        
        fiberstat = (petal['COADD_FIBERSTATUS'] == 0) 

        petals.append(petal[fiberstat & ~badfiber])
    
    return atable.vstack(petals) 


def rr_deep(tileid, release='everest'): 
    ''' redrock redshift success for given tile ID for cumulative coadds. 
    '''
    if release == 'everest': 
        dat_dir = '/global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/'
    else: 
        raise ValueError
    
    assert len(glob.glob(os.path.join(dat_dir, str(tileid), '*'))) == 1
    exp_dir = glob.glob(os.path.join(dat_dir, str(tileid), '*'))[0]
    date    = os.path.basename(exp_dir)

    fRR = lambda i: os.path.join(exp_dir, 'redrock-%i-%i-thru%s.fits' % (i, tileid, date))
    
    petals = [] 
    
    for i in range(10): 
        frr = fRR(i)

        if not os.path.isfile(frr): 
            print("  %s does not exist" % frr) 
            continue 

        hdu     = fits.open(frr)
        zbest   = atable.Table(hdu['REDSHIFTS'].data) # redrock output 
        fmap    = atable.Table(hdu['FIBERMAP'].data)  # fibermap data 
        tsnr2   = atable.Table(hdu['TSNR2'].data)     # TSNR2 data 
        
        # zbest and fmap should be row-matched starting from Everest :) 
        petal   = atable.join(zbest, fmap,  keys='TARGETID', join_type='left')
        petal   = atable.join(petal, tsnr2, keys='TARGETID', join_type='left')
        
        # bad fibers
        badfiber  = (petal['ZWARN'] & 2**9) != 0 # no data 
        badfiber |= (petal['ZWARN'] & 2**11) != 0 # poor positioning
        badfiber |= (petal['ZWARN'] & 2**7)  != 0 # broken
        
        fiberstat = (petal['COADD_FIBERSTATUS'] == 0) 

        petals.append(petal[fiberstat & ~badfiber])
    
    return atable.vstack(petals) 

def rr_deep_hp(tileid, release='everest', survey='sv1'): 
    ''' redrock redshift success for given tile ID for cumulative coadds. 
    '''
    
    assert survey in ['sv3', 'main']
    
    if release == 'everest': 
        # /sv3/bright/
        dat_dir = '/global/cfs/cdirs/desi/spectro/redux/everest/healpix/'
    else: 
        raise ValueError

    ff = open('/global/cfs/cdirs/desi/spectro/redux/everest/healpix/tilepix.json')
    tpix = json.load(ff)    
    ff.close()
    
    hps  = [tpix['{}'.format(tileid)]['{}'.format(petal)] for petal in range(10)]
    
    # flatten.
    hps  = np.array([item for sublist in hps for item in sublist])

    result = []
    
    for hp in hps.astype(str): 
        # e.g. /global/cfs/cdirs/desi/spectro/redux/everest/healpix/sv3/bright/260/26050/redrock-sv3-bright-26050.fits
        frr = dat_dir + '/{}/bright/{}/{}/redrock-{}-bright-{}.fits'.format(survey, hp[:-2], hp, survey, hp)

        if not os.path.isfile(frr): 
            print("  %s does not exist" % frr) 
            continue 

        hdu     = fits.open(frr)
        zbest   = atable.Table(hdu['REDSHIFTS'].data) # redrock output         
        fmap    = atable.Table(hdu['FIBERMAP'].data)  # fibermap data 
        tsnr2   = atable.Table(hdu['TSNR2'].data)     # TSNR2 data 
        
        # zbest and fmap should be row-matched starting from Everest :) 
        hp      = atable.join(zbest, fmap,  keys='TARGETID', join_type='left')
        hp      = atable.join(hp,    tsnr2, keys='TARGETID', join_type='left')
        
        # bad fibers
        badfiber  = (hp['ZWARN'] & 2**9) != 0 # no data 
        badfiber |= (hp['ZWARN'] & 2**11) != 0 # poor positioning
        badfiber |= (hp['ZWARN'] & 2**7)  != 0 # broken
        
        fiberstat = (hp['COADD_FIBERSTATUS'] == 0) 

        result.append(hp[fiberstat & ~badfiber])
    
    return atable.vstack(result) 

def set_zbest_exp_zsuccess(zbest_exp, exp_dX2=40.):
    zbest_exp  = atable.Table(zbest_exp, copy=True) 
    
    crit_zwarn = (zbest_exp['ZWARN'] == 0)
    crit_dchi2 = (zbest_exp['DELTACHI2']  > exp_dX2) 
    crit_stype = (zbest_exp['SPECTYPE'] != "STAR") & (zbest_exp['SPECTYPE'] != "QSO")# only galaxy spectra
    crit_z_lim = (zbest_exp['Z'] > 0.0) & (zbest_exp['Z'] < 0.6) # rough BGS redshift limit
    crit_z_err = (zbest_exp['ZERR'] < (0.0005 * (1. + zbest_exp['Z'])))
    
    dz_1pz     = np.abs(zbest_exp['Z_TRUE'] - zbest_exp['Z']) / (1. + zbest_exp['Z_TRUE'])
    crit_ztrue = (dz_1pz < 0.003)

    # combine all criteria
    zsuccess   = (crit_zwarn & crit_dchi2 & crit_stype & crit_z_lim & crit_z_err
                  & crit_ztrue & zbest_exp['DEEP_TRUE'])
    
    zbest_exp['ZSUCCESS']       = zsuccess
    zbest_exp['RR_ZWARN']       = zbest_exp['ZWARN']
    zbest_exp['RR_DELTACHI2']   = zbest_exp['DELTACHI2']
    zbest_exp['RR_SPECTYPE']    = zbest_exp['SPECTYPE']
    zbest_exp['RR_ZERR']        = zbest_exp['ZERR']
    zbest_exp['RR_Z']           = zbest_exp['Z']
    zbest_exp['RR_Z_DEEP']      = zbest_exp['Z_TRUE']
    
    return zbest_exp

def get_zbest_exp(tileid, expid, release='everest', survey='sv1', ext_cols=None, exp_dX2=40.): 
    ''' get redshift success rate for given exposure 
    '''
    # get redrock file for deep, which will be used as the truth table.
    _zbest_deep = rr_deep(tileid, release=release)

    if survey in ['sv3', 'main']:
        # use deepest coadds - healpixel (coadded across tiles).
        dp_tids     = _zbest_deep['TARGETID'].data
        
        # targets in hp overlapping tile.
        _zbest_deep = rr_deep_hp(tileid, release=release, survey=survey)
        
        # remove targets not in tile. 
        _zbest_deep = _zbest_deep[np.isin(_zbest_deep['TARGETID'], dp_tids)]
            
    crit_stype  = (_zbest_deep['SPECTYPE'] != "STAR") & (_zbest_deep['SPECTYPE'] != "QSO") # only galaxy spectra
    crit_z_lim  = (_zbest_deep['Z'] > 0.0) & (_zbest_deep['Z'] < 0.6) # rough BGS redshift limit
    crit_zwarn  = (_zbest_deep['ZWARN'] == 0)
    crit_dchi2  = (_zbest_deep['DELTACHI2']  > 40.) 
    crit_z_err  = (_zbest_deep['ZERR'] < (0.0005 * (1. + _zbest_deep['Z'])))

    cols        = ['TARGETID', 'Z', 'COADD_NUMEXP', 'COADD_EXPTIME', 'TSNR2_BGS']

    if ext_cols != None:
        cols += ext_cols

    # unique
    cols        = list(set(cols))

    # exclude spectra that are stars or outside of the BGS redshift range in
    # the redhsift success calculations. 
    keep        = crit_stype & crit_z_lim
    zbest_deep  = _zbest_deep[keep][cols]

    # redshift success of deep exposure --- this is the best we can possibly do
    deep_true   = (crit_zwarn & crit_dchi2 & crit_z_err)[keep]
    zbest_deep['DEEP_TRUE'] = deep_true

    # redrock redshifts from DEEP exposure will be used as true redshifts
    zbest_deep.rename_column('Z', 'Z_TRUE')
    zbest_deep.rename_column('COADD_NUMEXP',  'DEEP_NUMEXP')
    zbest_deep.rename_column('COADD_EXPTIME', 'DEEP_EXPTIME')
    zbest_deep.rename_column('TSNR2_BGS',     'DEEP_TSNR2_BGS')
    
    # get redrock output for SINGLE exposure for BGS_ANY targets
    zbest_exp = rr_exposure(tileid, expid, release=release)

    # only keep targets that are in the deep
    zbest_exp = atable.join(zbest_deep, zbest_exp, keys='TARGETID', join_type='inner')
    
    zbest_exp = set_zbest_exp_zsuccess(zbest_exp, exp_dX2=exp_dX2)
        
    # mags. 
    zbest_exp = zbest_exp[(zbest_exp['PHOTSYS'] != '') & (zbest_exp['PHOTSYS'] != 'G')]

    for band in ['g', 'r', 'z', 'W1', 'W2']:
        trans = mwdust_transmission(zbest_exp['EBV'].data, band, zbest_exp['PHOTSYS'].data.astype(str), match_legacy_surveys=False)

        zbest_exp['{}MAG_DRED'.format(band.upper())] = 22.5 - 2.5 * np.log10((zbest_exp['FLUX_{}'.format(band.upper())] / trans).clip(1e-16))

        if band == 'r':
            zbest_exp['FIBER_RMAG_DRED'] = 22.5 - 2.5 * np.log10((zbest_exp['FIBERFLUX_{}'.format(band.upper())] / trans).clip(1e-16))

    zbest_exp['FAINT_FIBCOL'] = (zbest_exp['ZMAG_DRED'] - zbest_exp['W1MAG_DRED']) - 3. / 2.5 * (zbest_exp['GMAG_DRED'] - zbest_exp['RMAG_DRED']) + 1.2
            
    return zbest_exp


def zsuccess_rate(prop, zsuccess_cond, range=None, nbins=20, bin_min=2):
    ''' measure the redshift success rate along with property `prop`

    :params prop: 
        array of properties (i.e. Legacy r-band magnitude) 

    :params zsuccess_cond:
        boolean array indicating redshift success 

    :params range: (default: None) 
        range of the `prop` 

    :params nbins: (default: 20) 
        number of bins to divide `prop` by 
    
    :params bin_min: (default: 2)  
        minimum number of objects in bin to exlcude it 

    :return wmean: 
        weighted mean of `prop` in the bins 

    :return e1: 
        redshift success rate in the bins

    :return ee1: 
        simple poisson error on the success rate
    '''
    h0, bins = np.histogram(prop, bins=nbins, range=range)
    hv, _ = np.histogram(prop, bins=bins, weights=prop)
    h1, _ = np.histogram(prop[zsuccess_cond], bins=bins)
    
    good = h0 > bin_min
    hv = hv[good]
    h0 = h0[good]
    h1 = h1[good]

    wmean = hv / h0 # weighted mean 
    rate = h1.astype("float") / (h0.astype('float') + (h0==0))
    e_rate = np.sqrt(rate * (1 - rate)) / np.sqrt(h0.astype('float') + (h0 == 0))
    return wmean, rate, e_rate
