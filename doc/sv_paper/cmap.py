import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def get_cmap(): 
    clrs = [] 

    for i in range(10): 
        clrs.append((0, 0, 0.1 * i, 1))
        #clrs.append((0.05 * i, 0.05 * i, 0.05 * i, 1))

    clrs = clrs[::-1]
    cmap = mpl.colors.ListedColormap(clrs)#[(0,0,0,0.), (0,0,0,0.95)])

    # sample the colormaps that you want to use. Use 128 from each so we get 256
    # colors in total
    clrs1 = plt.cm.Blues(np.linspace(0.5, 1, 5))
    clrs2 = plt.cm.Spectral_r(np.linspace(0.2, 1, 95))
    #clrs2 = plt.cm.hot_r(np.linspace(0.2, 0.8, 95))

    # combine them and build a new colormap
    colors = np.vstack((clrs2, clrs1))
    mymap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)
    
    return mymap

cmap = get_cmap()