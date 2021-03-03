import numpy as np


def is_badz(cat, dX2_lim=40., verbose=False, summary=False, cuts_list=None):
    
    badz     = np.zeros(len(cat)).astype(bool) 

#     labels   = ['ZWARN', 'DELTACHI2', 'SPECTYPE', 'ZRANGE', 'ZERR']
#     #labels   = ['ZWARN', 'ZERR']

#     cuts     = [cat['ZWARN'] > 0,\
#                 cat['DELTACHI2'] < dX2_lim,\
#                 cat['SPECTYPE'] == 'STAR',\
#                 (cat['Z'] < 0.0) | (cat['Z'] > 0.6),\
#                 cat['ZERR'] > (0.0005 * (1. + cat['Z']))
#                ]
    
    cuts     = {}
    cuts['ZWARN'] = cat['ZWARN'] > 0
    cuts['DELTACHI2'] = cat['DELTACHI2'] < dX2_lim
    cuts['SPECTYPE'] = cat['SPECTYPE'] == 'STAR'
    cuts['ZRANGE'] = (cat['Z'] < 0.0) | (cat['Z'] > 0.6)
    cuts['ZERR'] = cat['ZERR'] > (0.0005 * (1. + cat['Z']))

    isummary = {}
    
    if cuts_list is None:
        print('is None')
        #for label, cut in zip(labels, cuts):
        for label, cut in cuts.items():
            badz = badz | cut

            if verbose:
                print('LOST {} TO {} CUT.'.format(np.count_nonzero(cut), label))

            isummary[label] = np.count_nonzero(cut)
    else:
        print('is Not None')
        for label in cuts_list:
            cut = cuts[label]
            badz = badz | cut

            if verbose:
                print('LOST {} TO {} CUT.'.format(np.count_nonzero(cut), label))

            isummary[label] = np.count_nonzero(cut)

    if summary:
        return  badz, isummary, cuts

    else:
        return  badz, cuts