#!/usr/bin/env python
import glob
import numpy as np
import astropy.io.fits as fits

from   astropy.table import Table
from   redrock.rebin import trapz_rebin


todo      = glob.glob('spectra/lynx/coadd*')

for x in todo:
    coadd = fits.open(x)
    zbest = Table(fits.open('lynx/{}'.format(x.split('/')[-1].replace('coadd', 'zbest')))['ZBEST'].data)

    # coadd-6-80613-all.fits
    name  = x.split('/')[-1]
    parts = name.split('-')

    petal = parts[1]
    tile  = parts[2]
        
    wave  = coadd['BRZ_WAVELENGTH'].data
    flux  = coadd['BRZ_FLUX'].data
    tid   = coadd['FIBERMAP'].data['TARGETID']
    fstat = coadd['FIBERMAP'].data['FIBERSTATUS']
    
    # Sort by TARGETID order. 
    idx   = np.argsort(tid)
    flux  = flux[idx,:]

    # 
    zwarn = zbest['ZWARN'].data
    dchi2 = zbest['DELTACHI2'].data

    keep  = zbest['TARGETID'].data 

    assert  np.all(np.isin(keep, tid))     # Complete in targetid. 
    assert  np.all(keep == np.sort(keep))  # Targetid order.
    
    keep  = (zwarn == 0) & (dchi2 > 40.)

    print('\n\n')
    print('Petal {} on Tile {}'.format(petal, tile))
    print('{} failed fibers'.format(np.count_nonzero(fstat)))
    print('{} redshifts.'.format(len(zbest)))
    print('{} successful.'.format(np.count_nonzero(keep)))

    flux  = flux[keep,:]

    # Wave [A] and Flux 1.e-17 erg/s/cm2/A.
    result = np.c_[wave, flux.T]

    np.savetxt('spectra/qspec_in/bgs-{}-{}.txt'.format(petal, tile), result, header='# Wave [A]    FLUX [1.e-17 erg/s/cm2/A]')

print('\n\nDone.\n\n')
