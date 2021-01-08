#!/usr/bin/env python
from setuptools import find_packages, setup

__version__ = '0.0'

if __name__ == "__main__":
    setup(
        name='bgs_sv',
        version = __version__,
        description='python package for DESI BGS CMX + SV efforts',
        packages=['bgs_sv'],
        package_dir={'': 'py'}, 
        install_requires = ['numpy', 'matplotlib', 'scipy', 'fitsio'],
        include_package_data=True,
        package_data={'bgs_sv': ['dat/*.fits', 'dat/*.txt']},
        zip_safe=True,
    )
