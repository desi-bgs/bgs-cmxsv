#!/usr/bin/env python 
import os
import glob

'''
Coadd deep field cframe across cameras for input to e.g. desisim/scripts/quickspectra or specsim. 
'''

# Point to deep field coadds.
todo = [x.split('/')[-1] for x in glob.glob('lynx/coadd*.fits')]

for x in todo:
    cmd = 'desi_coadd_spectra -i lynx/{} --coadd-cameras -o spectra/lynx/{}'.format(x, x)

    # desi_coadd_spectra -i lynx/coadd-0-80613-all.fits --coadd-cameras -o spectra/lynx/coadd-0-80613-all.fits 
    print(cmd)

    os.system(cmd)

