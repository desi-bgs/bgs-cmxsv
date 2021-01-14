#!/usr/bin/env python                                                                                                                                                                                                                      
import os
import glob
import numpy as np
import astropy.io.fits as fits

from   astropy.table import Table


exps  = Table.read('bgs-cmxsv/py/bgs-cmxsv/dat/sv1-exposures.fits')

brite =  exps[ exps['GFA_SKY_MAG_AB_MED']   < 20.07]
brite = brite[brite['GFA_TRANSPARENCY_MED'] >  0.95]
brite = brite[brite['GFA_MOON_ZD_DEG_MIN']  <   90.]
brite = brite[brite['GFA_FWHM_ASEC_MED']    <   1.8]

print(brite)

# /global/homes/m/mjwilson/desi/SV1/spec-lrg-z0.7-zmag20.00.dat
todo         = glob.glob('spectra/qspec_in/*.txt')

np.random.seed(314)

for cond in brite:
    expid    = cond['EXPID']

    seeing   = cond['GFA_FWHM_ASEC_MED']
    airmass  = cond['GFA_AIRMASS_MED']
    exptime  = cond['EXPTIME']
    moonfrac = cond['GFA_MOON_ILLUMINATION_MED']
    moonalt  = 90. - cond['GFA_MOON_ZD_DEG_MED']
    moonsep  = cond['GFA_MOON_SEP_DEG_MED']
    
    x        = np.random.choice(todo, replace=True)
    
    out      = x.replace('in', 'out/{:08d}/'.format(expid)).replace('bgs', 'bgs-{:08d}'.format(expid)).replace('txt', 'fits')
    cmd      = 'quickspectra --program BRIGHT --seeing {:.2f} --airmass {:.2f} --exptime {:.1f} --moonfrac {:.3f} --moonalt {:.1f} --moonsep {:.1f} --source-type bgs\
                             --config_file /global/homes/m/mjwilson/desi/SV1/specsim/specsim/data/config/desi.yaml --i {} --o {}'.format(seeing, airmass, exptime, moonfrac, moonalt, moonsep, x, out)

    os.system(cmd)
        
# rm spectra/zbest-bgs.fits
# rm spectra/redrock-bgs.h5

# rrdesi /global/cscratch1/sd/mjwilson/desi/SV1/spectra/bgs.fits -o /global/cscratch1/sd/mjwilson/desi/SV1/spectra/redrock-bgs.h5  -z /global/cscratch1/sd/mjwilson/desi/SV1/spectra/zbest-bgs.fits
