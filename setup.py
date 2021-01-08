#!/usr/bin/env python
from setuptools import find_packages, setup

__version__ = '0.0'

PACKAGES = find_packages(where="py")

if __name__ == "__main__":
    setup(
        name='bgs-cmxsv',
        version = __version__,
        description='python package for DESI BGS CMX + SV efforts',
        packages=PACKAGES,
        package_dir={"": "py"},
        install_requires = ['numpy', 'matplotlib', 'scipy'],
        zip_safe=True,
    )
