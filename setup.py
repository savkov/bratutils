#!/usr/bin/env python

import os

from setuptools import setup, find_packages

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open('VERSION') as f:
    VERSION = f.read().splitlines()[0]

setup(
    name='bratutils',
    version=VERSION,
    packages=find_packages('src', exclude=('tests',)),
    package_dir={'': 'src'},
    zip_safe=True,
    include_package_data=False,
    description='brat utilities',
    author='Sasho Savkov',
    license='GNU',
    long_description=(
        'https://github.com/savkov/bratutils'
    )
)
