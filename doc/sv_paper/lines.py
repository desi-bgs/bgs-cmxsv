import astropy.table as atable

def get_lines(survey='sv3', release='everest'):
    # FASTPHOT, METADATA. 
    fastphot = atable.Table.read('/global/cfs/cdirs/desi/spectro/fastspecfit/{}/catalogs/fastphot-{}-{}-bright.fits'.format(release, release, survey), 
                                 'FASTPHOT')
    fp_meta  = atable.Table.read('/global/cfs/cdirs/desi/spectro/fastspecfit/{}/catalogs/fastphot-{}-{}-bright.fits'.format(release, release, survey), 
                                 'METADATA')

    # FASTSPEC, METADATA.
    fastspec = atable.Table.read('/global/cfs/cdirs/desi/spectro/fastspecfit/{}/catalogs/fastspec-{}-{}-bright.fits'.format(release, release, survey), 
                                 'FASTSPEC')
    fs_meta  = atable.Table.read('/global/cfs/cdirs/desi/spectro/fastspecfit/{}/catalogs/fastspec-{}-{}-bright.fits'.format(release, release, survey), 
                                 'METADATA')
    
    return fastphot, fp_meta, fastspec, fs_meta

def lines_match(sv_gals, release='everest', survey='sv3'):
    fastphot, fp_meta, fastspec, fs_meta = get_lines(release=release, survey=survey)

    sv_gals   = atable.Table(sv_gals, copy=True)
    
    linelist  = ['OII_3726', 'OII_3729', 'HBETA', 'OIII_4363', 'OIII_4959', 'HALPHA', 'SII_6731', 'SII_6716']
    linelist  = ['{}_FLUX'.format(line) for line in linelist]
    linelist += ['TARGETID']
    
    return  atable.join(sv_gals, fastspec[linelist], keys=['TARGETID'], join_type='left')
    