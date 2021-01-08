# bgs-cmxsv
BGS survey planning during Commissioning and Survey Validation phases of the
DESI survey.

## installation
If you're ssh-ing into NERSC 
```
git clone https://github.com/desi-bgs/bgs-cmxsv.git
cd bgs-cmxsv
pip install . 
```

If you're on NERSC jupyter notebook run the following in the notebook 
```
!python -m pip install git+https://github.com/desi-bgs/bgs-cmxsv.git
```

## usage
```
from bgs_sv import sv1 # module for reading in sv1 data 
```
