import astropy.table as atable

def get_gama():
    lgama = atable.Table.read('/global/cfs/cdirs/desi/target/analysis/truth/dr9.0/south/matched/ls-dr9.0-GAMA-DR3-SpecObj-match.fits')
    zgama = atable.Table.read('/global/cfs/cdirs/desi/target/analysis/truth/dr9.0/south/matched/GAMA-DR3-SpecObj-match.fits')
    
    lgama['GAMA_NAME'] = zgama['GAMA_NAME'].data
    lgama['GAMA_SPECID']    = zgama['SPECID'].data
    lgama['GAMA_SURVEY']    = zgama['SURVEY'].data
    lgama['GAMA_Z']     = zgama['Z'].data
    lgama['GAMA_NQ']        = zgama['NQ'].data

    lgama.rename_column('OBJID', 'BRICK_OBJID')
    
    return lgama['GAMA_NAME', 'GAMA_SPECID', 'GAMA_SURVEY', 'GAMA_Z', 'GAMA_NQ', 'BRICKID', 'BRICK_OBJID']
    
def gama_match(sv_gals):
    lgama   = get_gama()
    sv_gals = atable.Table(sv_gals, copy=True)
        
    return atable.join(sv_gals, lgama, keys=['BRICK_OBJID', 'BRICKID'], join_type='left')