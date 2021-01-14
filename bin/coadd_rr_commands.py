# Get the a list of science exposures for a specficic date and tile
from    __future__ import division, print_function

import  ephem
import  fitsio
import  healpy
import  desisurvey.config
import  sys, os, glob, time, warnings
import  numpy                as     np
import  astropy.units        as     u
import  desisurvey.utils     as     dutils

from    pathlib              import Path
from    astropy.time         import Time
from    astropy.table        import Table, vstack, hstack
from    astropy.coordinates  import Angle
from    desisurvey.utils     import get_airmass
from    desiutil             import dust
from    astropy.coordinates  import SkyCoord

# os.system('source ./env.sh')

date          = '20210113'
redux_dir     = '/global/cfs/cdirs/desi/spectro/redux/blanc/'
output_dir    = '/global/homes/m/mjwilson/desi/SV1/spectra/exposures/'

# number of exposures in a coadded; 1 for single-exposure coadd                                                                                                                                                                             
ALL           = False   # Overrides n_exp.

n_exp         = 1
n_node        = 4
nside         = 512

overwrite     = True
archetypes    = False
verbose       = False

#
cond          = Table.read('bgs-cmxsv/py/bgs-cmxsv/dat/sv1-exposures.fits')
cond          = cond[cond['TARGETS'] == 'BGS+MWS']

nights        = np.unique(cond['NIGHT'].data).astype(str)
tiles         = np.unique(cond['TILEID'].data)

petals        = list(range(10))

print('\n\n\t---  TILES  ---')

for tile in tiles:
  print('\t    {}'.format(tile))

print('\n\t---  NIGHTS  ---')

for night in nights:
  print('\t    {}'.format(night))

print('\n')

if not os.path.isfile('./spectra/bgs_allcframes_{}.fits'.format(date)):
    ################################## Get list of exposures ##################################
    exposure_dir_list    = []
    
    for obsdate in nights:
        exposure_dir_list += glob.glob(os.path.join(redux_dir, 'exposures', obsdate, '*'))

    ################################## Get list of all science cframes ##################################  
    cframe_list = []

    # Get a list of all science exposures.
    for exposure_dir in exposure_dir_list:
        cframe_list_tmp = glob.glob(os.path.join(exposure_dir, 'cframe-*'))

        if len(cframe_list_tmp) > 0:
            if tiles is None:
                cframe_list += cframe_list_tmp

        else:
            # only need to check one cframe file in the exposure
            with fitsio.FITS(cframe_list_tmp[0]) as f:
                if f[0].read_header()['TILEID'] in tiles:
                    cframe_list += cframe_list_tmp
                    
    cframe_list           = sorted(cframe_list)

    # Gather exposure/petal information
    cframes               = Table()

    cframes['cframe']     = np.array(cframe_list)
    cframes['night']      = '                                        '
    cframes['mjd']        = np.zeros(len(cframes), dtype=np.float) - 1.0 
    cframes['lat']        = np.zeros(len(cframes), dtype=np.float)
    cframes['lon']        = np.zeros(len(cframes), dtype=np.float)
    cframes['elv']        = np.zeros(len(cframes), dtype=np.float)
    cframes['tileid']     = np.zeros(len(cframes), dtype=int)
    cframes['expid']      = np.zeros(len(cframes), dtype=int)
    cframes['exptime']    = np.zeros(len(cframes), dtype=np.float)
    cframes['camera']     = '                                        '
    cframes['program']    = '                                        '
    cframes['petal_loc']  = -1 * np.ones(len(cframes), dtype=np.int32)
    cframes['specgrph']   = np.zeros(len(cframes), dtype=np.int)
    cframes['ra']         = np.zeros(len(cframes), dtype=np.float)
    cframes['dec']        = np.zeros(len(cframes), dtype=np.float)

    for index, cframe in enumerate(cframes['cframe']):
        with fitsio.FITS(cframe) as f:
            header                         = f[0].read_header()

            if verbose:
                print(header)
           
            cframes['mjd'][index]          = header['MJD-OBS']
            cframes['night'][index]        = header['NIGHT']
            cframes['tileid'][index]       = header['TILEID']
            cframes['expid'][index]        = header['EXPID']
            cframes['camera'][index]       = header['CAMERA'].strip()[0]
            cframes['petal_loc'][index]    = int(header['CAMERA'].strip()[1])
            cframes['program'][index]      = header['PROGRAM']
            cframes['lat'][index]          = header['OBS-LAT']
            cframes['lon'][index] 	       = header['OBS-LONG']
            cframes['elv'][index] 	       = header['OBS-ELEV']
            cframes['exptime'][index]      = header['EXPTIME']
            cframes['ra'][index]           = header['SKYRA']
            cframes['dec'][index]          = header['SKYDEC']
            cframes['specgrph'][index]     = header['SPECGRPH']

    # Sanity check: each petal must have three cframe files.
    for expid in np.unique(cframes['expid']):
        mask_expid = cframes['expid'] == expid
    
        for petal_loc in petals:
            mask = mask_expid & (cframes['petal_loc'] == petal_loc)

            if (np.sum(mask) > 0) & (np.sum(mask) != 3):
                raise  ValueError('EXPID {} PETAL_LOC {} has only {} cframes files'.format(expid, petal_loc, np.sum(mask)))

    print('\n\n')
    
    uids, cnts = np.unique(cframes['expid'], return_counts=True)

    cframes.sort(('tileid', 'petal_loc'))

    ##  cframes.pprint(max_width=-1)

    cframes.write('./spectra/bgs_allcframes_{}.fits'.format(date), format='fits', overwrite=True)

## 
cframes = Table.read('./spectra/bgs_allcframes_{}.fits'.format(date))

## 
output_argument_list = []

if not ALL:
  if (not overwrite) and (os.path.isfile(output_dir + "/scripts/commands_coadd_nexp_{}_{}.sh".format(n_exp, night)) | os.path.isfile(output_dir + "/scripts/commands_rr_nexp_{}_{}.sh".format(n_exp, night))):
    raise  ValueError('Overwrite=True required to remove exisiting files.')

  output_file   = open(output_dir + "/scripts/commands_coadd_nexp_{}_{}.sh".format(n_exp, night), "w")
  output_rrfile = open(output_dir + "/scripts/commands_rr_nexp_{}_{}.sh".format(n_exp, night),    "w")
  
else:
  if (not overwrite) and (os.path.isfile(output_dir + "/scripts/commands_coadd_allexp_{}.sh".format(night)) | os.path.isfile(output_dir + "/scripts/commands_rr_allexp_{}.sh".format(night))):
    raise  ValueError('Overwrite=True required to remove exisiting files.')
  
  output_file   = open(output_dir + "/scripts/commands_coadd_allexp_{}.sh".format(night), "w")
  output_rrfile = open(output_dir + "/scripts/commands_rr_allexp_{}.sh".format(night),    "w")

Path(output_dir).mkdir(parents=True, exist_ok=True)
  
##
for night in nights:
  nightframes   = cframes[cframes['night'] == night]

  for tileid in np.unique(nightframes['tileid']):    
    for petal_loc in petals: 
        mask                    = (cframes['tileid'] == tileid) & (cframes['petal_loc'] == petal_loc) & (cframes['night'] == night)

        # Check b camera for simplicity, only reduce if all three cameras are present.  
        mask                   &= (cframes['camera'] == 'b')

        cframe1                 = cframes[mask]

        if not ALL:
          if (np.count_nonzero(mask) < n_exp):
            print('\t\t {} exposures is not enough for TILEID {}, NIGHT {} PETAL_LOC {} for a {} exp. coadd.'.format(np.sum(mask), tileid, night, petal_loc, n_exp))
            continue

          else:
            print('\t {} exposures for TILEID {}, NIGHT {} PETAL_LOC {} (for a {} exp. coadd).'.format(np.sum(mask), tileid, night, petal_loc, n_exp))
          
          # Skip the exposures that do not make the split.
          cframe1               = cframe1[:len(cframe1) - len(cframe1) % (n_exp)]

          nsplit                = len(cframe1)//(n_exp)

          subset_split          = np.split(np.arange(len(cframe1)), nsplit)

        else:
          if (np.sum(mask) == 0):
            print('\n# No exposures for TILEID {}, PETAL_LOC {}.\n'.format(tileid, petal_loc))
            continue
          
          nsplit                = 1
          subset_split          = np.split(np.arange(len(cframe1)), nsplit)
          
        for subset_index in range(len(subset_split)):
            subset              = cframe1[subset_split[subset_index]]
            input_argument      = ''

            for index in range(len(subset)):
                exposure_dir    = os.path.dirname(subset['cframe'][index])

                # .replace(redux_dir, '$REDUXDIR/')
                input_argument += os.path.join(exposure_dir, 'cframe-[brz]{}-*.fits ').format(petal_loc)

            if (not ALL) & (n_exp == 1):
                exposure        = os.path.basename(exposure_dir)

                output_argument = os.path.join(output_dir, 'NEXP{}'.format(n_exp), str(tileid), night, 'coadd-{}-{}-{}.fits'.format(night, petal_loc, exposure))

            elif not ALL:
                output_argument = os.path.join(output_dir, 'NEXP{}'.format(n_exp), str(tileid), night, 'coadd-{}-{}-{}exp-subset-{}.fits'.format(night, petal_loc, n_exp, subset_index))

            else:
                output_argument = os.path.join(output_dir, 'ALL', str(tileid), night, 'coadd-{}-{}-allexp.fits'.format(night, petal_loc))
                
            output_argument_list.append(output_argument)

            if os.path.isfile(output_argument) and (not overwrite):
                print('\nWarninig: {} already exists!\n'.format(output_argument))
                continue

            # --coadd-cameras  
            output_file.write('time desi_coadd_spectra -i {} -o {}\n'.format(input_argument, output_argument))

##
output_file.close()
        
for output_argument in output_argument_list:
    rrdesi_argument_redrock = output_argument.replace('coadd', 'redrock').replace('.fits', '.h5')
    rrdesi_argument_zbest   = output_argument.replace('coadd', 'zbest')

    cmd                     = 'srun -N {} -n {} -c 2 rrdesi_mpi -o {} -z {} {}'.format(n_node, 32 * n_node, rrdesi_argument_redrock, rrdesi_argument_zbest, output_argument)

    if archetypes:
      cmd                   = cmd + ' --archetypes /global/common/software/desi/cori/desiconda/20190804-1.3.0-spec/code/redrock-archetypes/master/\n'    

    else:
      cmd                  += '\n'
      
    output_rrfile.write(cmd)

print('\n\nWritten commands bash script to {}/scripts.'.format(output_dir))
    
print('\n\nDone.\n\n')
