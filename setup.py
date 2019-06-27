#!/usr/bin/env python

import os

from distutils.core import setup
from setuptools import find_packages

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='bratutils',
    version='0.3.0.dev0',
    packages=find_packages('src', exclude=('tests',)),
    url='https://github.com/savkov/transcriptor',
    author='Sasho Savkov',
    author_email='me@sasho.io',
    package_dir={'': 'src'},
    zip_safe=True,
    include_package_data=False,
    description='brat utilities',
    license='MIT',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Operating System :: Unix',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering'
    ]
)
