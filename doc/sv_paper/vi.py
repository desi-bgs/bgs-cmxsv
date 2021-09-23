import numpy as np
import astropy.table as atable

def get_vi():
    cols = ['TARGETID', 'best_z', 'best_quality', 'best_spectype']
    
    aa = atable.Table.read('/global/cfs/cdirs/desi/sv/vi/TruthTables/Andes/BGS/Truth_table_Andes_reinspection_BGS_66003_20200315_v1.csv')

    for x in ['z', 'quality', 'spectype']:
        aa.rename_column('best {}'.format(x, x), 'best_{}'.format(x, x))

    aa = aa[cols]
    
    bb = atable.Table.read('/global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/BGS/desi-vi_BGS_tile80613_nightdeep_merged_all_210202.csv')[cols]
    
    cc = atable.Table.read('/global/cfs/cdirs/desi/sv/vi/TruthTables/Cascades/BGS/desi-vi_SV_cascades_combination_BGS_all_210521.csv')[cols]
    
    vi = atable.vstack([aa, bb, cc])
    
    vi = vi[vi['best_quality'] >= 2.5]
    
    vi.sort('best_quality')
    
    # 4 targets with repeated VI. 
    vi = atable.unique(vi, keys='TARGETID', keep='last')
    
    vi.rename_column('best_z', 'VI_Z')
    vi.rename_column('best_quality', 'VI_Q')
    vi.rename_column('best_spectype', 'VI_SPECTYPE')
    
    return  vi

def vi_match(sv_gals):
    vi      = get_vi()
    sv_gals = atable.Table(sv_gals, copy=True)
        
    return atable.join(sv_gals, vi, keys=['TARGETID'], join_type='left')