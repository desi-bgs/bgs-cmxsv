'''

module for interfacing with SV data 

'''
import os 
import glob
import yaml
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

        # zbest and fmap should be row-matched starting from Everest :) 
        petal = atable.join(zbest, fmap, keys='TARGETID', join_type='left')
        
        # bad fibers
        badfiber = (petal['ZWARN'] & 2**9) != 0 # no data 
        fiberstat = (petal['COADD_FIBERSTATUS'] == 0) 

        petals.append(petal[fiberstat & ~badfiber])
    
    return atable.vstack(petals) 


def rr_deep(tileid, release='everest'): 
    ''' redrock redshift success for given tile ID for cumulative coadds 
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

        # zbest and fmap should be row-matched starting from Everest :) 
        petal = atable.join(zbest, fmap, keys='TARGETID', join_type='left')
        
        # bad fibers
        badfiber = (petal['ZWARN'] & 2**9) != 0 # no data 
        fiberstat = (petal['COADD_FIBERSTATUS'] == 0) 

        petals.append(petal[fiberstat & ~badfiber])
    
    return atable.vstack(petals) 


def get_exp_zsuccess(tileid, expid, release='everest'): 
    ''' get redshift success rate for given exposure 
    '''
    # get redrock file for deep, which will be used as the truth table.
    _zbest_deep = rr_deep(tileid, release=release)
    
    crit_zwarn = (_zbest_deep['ZWARN'] == 0)
    crit_dchi2 = (_zbest_deep['DELTACHI2']  > 40.) 
    crit_stype = (_zbest_deep['SPECTYPE'] != "STAR") # only galaxy spectra
    crit_z_lim = (_zbest_deep['Z'] > 0.0) & (_zbest_deep['Z'] < 0.6) # rough BGS redshift limit
    crit_z_err = (_zbest_deep['ZERR'] < (0.0005 * (1. + _zbest_deep['Z'])))
    
    truth_table = crit_dchi2 & crit_stype & crit_z_lim & crit_z_err
    zbest_deep = _zbest_deep[truth_table]['TARGETID', 'Z', 'COADD_NUMEXP', 'COADD_EXPTIME']
    # redrock redshifts from DEEP exposure will be used as true redshifts
    zbest_deep.rename_column('Z', 'Z_TRUE')
    zbest_deep.rename_column('COADD_NUMEXP', 'DEEP_NUMEXP')
    zbest_deep.rename_column('COADD_EXPTIME', 'DEEP_EXPTIME')
    
    # get redrock output for SINGLE exposure for BGS_ANY targets
    zbest_exp = rr_exposure(tileid, expid, release=release)

    # only keep targets that are in the deep
    zbest_exp = atable.join(zbest_deep, zbest_exp, keys='TARGETID', join_type='inner')
    
    crit_zwarn = (zbest_exp['ZWARN'] == 0)
    crit_dchi2 = (zbest_exp['DELTACHI2']  > 15.)#40.) 
    crit_stype = (zbest_exp['SPECTYPE'] != "STAR") # only galaxy spectra
    crit_z_lim = (zbest_exp['Z'] > 0.0) & (zbest_exp['Z'] < 0.6) # rough BGS redshift limit
    crit_z_err = (zbest_exp['ZERR'] < (0.0005 * (1. + zbest_exp['Z'])))
    
    dz_1pz = np.abs(zbest_exp['Z_TRUE'] - zbest_exp['Z']) / (1. + zbest_exp['Z_TRUE'])
    crit_ztrue = (dz_1pz < 0.003)

    # combine all criteria
    zsuccess = crit_zwarn & crit_dchi2 & crit_stype & crit_z_lim & crit_z_err & crit_ztrue
    zbest_exp['ZSUCCESS'] = zsuccess
    zbest_exp['RR_ZWARN']       = zbest_exp['ZWARN']
    zbest_exp['RR_DELTACHI2']   = zbest_exp['DELTACHI2']
    zbest_exp['RR_SPECTYPE']    = zbest_exp['SPECTYPE']
    zbest_exp['RR_ZERR']        = zbest_exp['ZERR']
    zbest_exp['RR_Z']           = zbest_exp['Z']
    zbest_exp['RR_Z_DEEP']      = zbest_exp['Z_TRUE']
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
