#!/usr/bin/env python
import os
import numpy as np
import astropy.io.fits as fits
import pylab as pl
import astropy

from   pathlib import Path
from   astropy.table import Table, vstack, join 

# https://github.com/desihub/desitarget/blob/master/py/desitarget/sv1/data/sv1_targetmask.yaml
from   desitarget.sv1.sv1_targetmask import desi_mask as sv1_desi_mask
from   desitarget.sv1.sv1_targetmask import bgs_mask as sv1_bgs_mask

# [BGS_FAINT,           0, "BGS faint targets",              {obsconditions: BRIGHT|GRAY|DARK}]
# [BGS_BRIGHT,          1, "BGS bright targets",             {obsconditions: BRIGHT}]
# [BGS_FAINT_EXT,       2, "BGS faint extended targets",     {obsconditions: BRIGHT}]
# [BGS_LOWQ,            3, "BGS low quality targets",        {obsconditions: BRIGHT}]
# [BGS_FIBMAG,          4, "BGS fiber magnitude targets",    {obsconditions: BRIGHT}]

# https://github.com/desihub/desitarget/blob/master/doc/nb/target-selection-bits-and-bitmasks.ipynb

exps  = Table.read('bgs-cmxsv/py/bgs-cmxsv/dat/sv1-exposures.fits')
exps  = exps[exps['TARGETS'] == 'BGS+MWS']

tiles = np.unique(exps['TILEID'].data)

# TILEID  NIGHT  EXPID
exps           = exps['TILEID', 'NIGHT', 'TILERA', 'TILEDEC', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']
exps['NEXP']   = 1

tinfo          = astropy.table.unique(exps['TILEID', 'TILERA', 'TILEDEC'], keys='TILEID')

exps_grouped   = exps.group_by(['TILEID', 'NIGHT'])
exps_binned    = exps_grouped.groups.aggregate(np.sum)
exps           = exps_binned

exps           = exps['TILEID', 'NIGHT', 'NEXP', 'EXPTIME', 'B_DEPTH', 'R_DEPTH', 'Z_DEPTH']
exps           = join(exps, tinfo, keys='TILEID', join_type='left')
  
exps.sort('TILEID')

exps['NBGSA'] = 0
exps['NBGSW'] = 0
exps['NBGSZ'] = 0

nights        = exps['NIGHT'].data.astype(np.int)
in_blanc      = nights <= 20201223

exps          = exps[in_blanc]
nights        = np.unique(exps['NIGHT'].data).astype(np.int)

# print(nights)

truth_cache = {}

def truth_test(tileid, cat):
  if tileid not in truth_cache.keys():  
    truth       = Table.read('/global/homes/m/mjwilson/desi/SV1/spectra/truth/bgs_deep_truth_{}.fits'.format(tileid))
    truth.sort('TARGETID')

    truth_cache[tileid] = truth

  else:
    truth    = truth_cache[tileid]

  in_cat     = np.isin(truth['TARGETID'], cat['TARGETID'])
  truth      = truth[in_cat]
    
  # Both sorted by targetid.
  for x in [truth, cat]:
    assert  np.all(np.argsort(x) == np.arange(len(x)))
    
  #
  in_truth   = np.isin(cat['TARGETID'], truth['TARGETID'])

  badz       = np.zeros(len(cat)).astype(bool)

  if (len(cat[in_truth]) > 0) & (len(truth) > 0):  
      badz[in_truth] = np.abs(cat['Z'][in_truth] - truth['Z']) > 2. * cat['ZERR'][in_truth]

  return  badz
  
  
for tileid in tiles:
  nights      = np.unique(exps[exps['TILEID'] == tileid]['NIGHT'])

  for night in nights:   
    for petal in range(10):
        # /global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/20201212/fba-080605.fits
        path         = '/global/homes/m/mjwilson/blanc/tiles/{}/{}/zbest-{}-{}-{}.fits'.format(tileid, night, petal, tileid, night)

        # print(path)
        
        if os.path.isfile(path):            
            infile   = fits.open(path)

            zbest    = infile['ZBEST'].data
            fmap     = infile['FIBERMAP'].data
            
            zbest    = Table(zbest)
            fmap     = Table(fmap)

            print('\n\n{} \t {} \t ({} \t {})'.format(tileid, petal, len(zbest), len(fmap)))
            
            tinfo    = fmap['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'FLUX_R', 'FIBERFLUX_R', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET'] 
            tinfo    = astropy.table.unique(tinfo, keys='TARGETID')
            
            # cols = fmap.dtype.names
            # cols = [x for x in cols if 'DESI_TARGET' in x]
            # print(tileid, cols)

            deep     = join(zbest, tinfo, keys='TARGETID', join_type='left')
            
            # FIBERSTATUS cut.                                                                                                                                                                                  
            # fibstat  = fmap['TARGETID', 'FIBERSTATUS']
            # fibstat  = fibstat.group_by('TARGETID')
            # fibstat  = fibstat.groups.aggregate(np.all)

            # badfiber = fibstat['FIBERSTATUS']

            # deep     = join(deep, fibstat, keys='TARGETID', join_type='left')

            assert  len(deep) == 500

            # Limit to BGS.                                                                                                                                                                                                       
            badbgs   = (deep['SV1_DESI_TARGET'] & sv1_desi_mask['BGS_ANY']) == 0

            nbright  = (deep['SV1_BGS_TARGET'] & sv1_bgs_mask['BGS_BRIGHT']) == 0
            nfaint   = (deep['SV1_BGS_TARGET'] & sv1_bgs_mask['BGS_FAINT'])  == 0
            
            # Limit to faint | bright bgs only.
            badbgs   = badbgs | (nbright & nfaint)
            
            # Limit to bright bgs only.                                                                                                                                                                                               
            badbgs   = badbgs | nbright
            
            #
            deep     = deep[~badbgs]
            
            # Assigned BGS, e.g. 3259 for 80614 here: https://data.desi.lbl.gov/desi/users/raichoor/fiberassign-sv1/sv1-per-tile/index.html#tile-nexp-design
            exps['NBGSA'][exps['TILEID'] == tileid] += len(deep)
            
            # https://github.com/desihub/redrock/blob/master/py/redrock/zwarning.py
            # NODATA for 'malfunctioning' positioner.
            deep['NODATA'] = (deep['ZWARN'] & 2**9) != 0 
            
            # print(np.count_nonzero(badfiber), np.count_nonzero(deep['NODATA']))

            # Limit to working fibers.
            badfiber = deep['NODATA']
            deep     = deep[~badfiber]

            exps['NBGSW'][exps['TILEID'] == tileid] += len(deep)
            
            # Good redshifts (in deep)
            badz     = (deep['ZWARN'] > 0) | (deep['DELTACHI2'] < 40.) | (deep['SPECTYPE'] == 'STAR')
            badz     = badz | (deep['Z'] < 0.0) | (deep['Z'] > 0.6)
            badz     = badz | (deep['ZERR'] > (0.0005 * (1. + deep['Z'])))
            
            # Test against deep 'truth'.
            # badz   = badz | truth_test(tileid, deep)

            #
            deep['BGS_SUCCESS'] = ~badz

            Path('spectra/bgs-zbest/{}/{}'.format(tileid, night)).mkdir(parents=True, exist_ok=True)
            
            deep.write('spectra/bgs-zbest/{}/{}/bgs-zbest-{}-{}-{}.fits'.format(tileid, night, petal, tileid, night))
            
            deep     = deep[~badz]
                        
            for label, lost in zip(['FIBER', 'BGS', 'Z'], [badfiber, badbgs, badz]):
                print('LOST {} on {} cut'.format(np.count_nonzero(lost), label))
                
            exps['NBGSZ'][exps['TILEID'] == tileid] += len(deep)

        else:
            print('Failed to retrieve {}.'.format(path)) 
  
exps['BGSSUCCESS_%'] = ['{:.2f}'.format(100. * x['NBGSZ'] / x['NBGSW']) for x in exps]

exps.sort('BGSSUCCESS_%')

print('\n\n')

exps.pprint(max_lines=-1)

print('\n\nDone.\n\n')
