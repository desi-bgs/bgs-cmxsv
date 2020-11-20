'''


module for reading in different BGS spectral simulations. Currently the module
includes: 
    * `simulated_GAMA_source_spectra`: read simulated GAMA-matched fiber-magnitude 
    scaled BGS source spectra


author: ChangHoon Hahn


'''
import os 
import h5py 
import numpy as np 

assert os.environ['NERSC_HOST'] == 'cori'


dir_spec_sim='/global/cfs/cdirs/desi/users/chahah/bgs_spec_sims' # directory with BGS spectral simulations 


def simulated_GAMA_source_spectra(emlines=True): 
    ''' read GAMA-matched fiber-magnitude scaled BGS source spectra 
    These source spectra are created for GAMA objects. their spectra is 
    constructed from continuum that's template matched to the broadband
    colors and emission lines from GAMA data (properly flux calibrated). 
    Then the spectra is scaled down to the r-band fiber magnitude. They 
    therefore do not require fiber acceptance fractions. 
    '''
    fsource = os.path.join(dir_spec_sim, 
            'GALeg.g15.sourceSpec%s.1000.seed0.hdf5' % ['.noemission', ''][emlines])

    if not os.path.isfile(fsource): 
        from feasibgs import catalogs as Cat 
        from feasibgs import forwardmodel as FM 
        seed = 0 
        np.random.seed(seed) 
        # read in GAMA-Legacy catalog with galaxies in both GAMA and Legacy surveys
        cata = Cat.GamaLegacy()
        gleg = cata.Read('g15', dr_gama=3, dr_legacy=7, silent=True)  
        
        # extract meta-data of galaxies 
        redshift        = gleg['gama-spec']['z']
        absmag_ugriz    = cata.AbsMag(gleg, kcorr=0.1, H0=70, Om0=0.3, galext=False) # ABSMAG k-correct to z=0.1
        r_mag_apflux    = 22.5 - 2.5 * np.log10(gleg['legacy-photo']['apflux_r'][:,1])) # aperture flux
        r_mag_gama      = gleg['gama-photo']['r_petro'] # r-band magnitude from GAMA (SDSS) photometry
        ha_gama         = gleg['gama-spec']['ha_flux'] # halpha line flux

        ngal = len(redshift) # number of galaxies
        vdisp = np.repeat(100.0, ngal) # velocity dispersions [km/s]

        # match GAMA galaxies to templates 
        bgs3 = FM.BGStree()
        match = bgs3._GamaLegacy(gleg)
        hasmatch = (match != -999)
        criterion = hasmatch 
        
        # randomly pick a few more than 5000 galaxies from the catalog that have 
        # matching templates because some of the galaxies will have issues where the 
        # emission line is brighter than the photometric magnitude.  
        subsamp = np.random.choice(np.arange(ngal)[criterion], int(1.1 * 1000), replace=False) 

        # generate noiseless spectra for these galaxies 
        s_bgs = FM.BGSsourceSpectra(wavemin=1500.0, wavemax=15000) 
        # emission line fluxes from GAMA data  
        if emlines: 
            emline_flux = s_bgs.EmissionLineFlux(gleg, index=subsamp, dr_gama=3, silent=True) # emission lines from GAMA 
            mag_em = r_mag_gama[subsamp]
        else: 
            emline_flux = None 
            mag_em = None 

        flux, wave, magnorm_flag = s_bgs.Spectra(
                r_mag_apflux[subsamp], 
                redshift[subsamp],
                vdisp[subsamp], 
                seed=1, 
                templateid=match[subsamp], 
                emflux=emline_flux, 
                mag_em=mag_em, 
                silent=True)

        # only keep 1000 galaxies
        isubsamp = np.random.choice(np.arange(len(subsamp))[magnorm_flag], 1000, replace=False) 
        subsamp = subsamp[isubsamp]
        
        # save to file  
        fsub = h5py.File(fsource, 'w') 
        fsub.create_dataset('zred', data=redshift[subsamp])
        fsub.create_dataset('absmag_ugriz', data=absmag_ugriz[:,subsamp]) 
        fsub.create_dataset('r_mag_apflux', data=r_mag_apflux[subsamp]) 
        fsub.create_dataset('r_mag_gama', data=r_mag_gama[subsamp]) 
        for grp in gleg.keys(): 
            group = fsub.create_group(grp) 
            for key in gleg[grp].keys(): 
                group.create_dataset(key, data=gleg[grp][key][subsamp])
        fsub.create_dataset('flux', data=flux[isubsamp, :])
        fsub.create_dataset('wave', data=wave)
        fsub.close()

    # read in source spectra
    source = h5py.File(fsource, 'r')

    wave_s = source['wave'][...]
    flux_s = source['flux'][...]
    
    meta = {} 
    for k in ['r_mag_apflux', 'r_mag_gama', 'zred', 'absmag_ugriz']: 
        meta[k] = source[k][...]
    meta['r_mag'] = 22.5 - 2.5 * np.log10(source['legacy-photo']['flux_r'][...])
    source.close()

    return wave_s, flux_s, meta
