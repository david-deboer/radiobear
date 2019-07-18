#! /usr/bin/env python
# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import absolute_import, division, print_function

import os
import glob
import io
from setuptools import setup, find_packages

from radiobear import version

with io.open('README.md', 'r', encoding='utf-8') as readme_file:
    readme = readme_file.read()

setup_args = {
    'name': "radiobear",
    'description': "radioBEAR:  radio version BErkeley Atmosphere Radiative transfer",
    'long_description': readme,
    'url': "https://github.com/david-deboer/radiobear",
    'license': "BSD",
    'author': "David DeBoer",
    'author_email': "ddeboer@berkeley.edu",
    'version': version.VERSION,
    'packages': ['radiobear', 'radiobear.plotting',
                 'radiobear.Jupiter', 'radiobear.Saturn', 'radiobear.Uranus', 'radiobear.Neptune',
                 'radiobear.Constituents',
                 'radiobear.Constituents.clouds', 'radiobear.Constituents.co', 'radiobear.Constituents.h2',
                 'radiobear.Constituents.h2o', 'radiobear.Constituents.h2s', 'radiobear.Constituents.nh3',
                 'radiobear.Constituents.ph3'],
    'scripts': glob.glob('scripts/*'),
    'include_package_data': True,
    'install_requires': ["six", "numpy", "matplotlib", "scipy"],
    'classifiers': ["Development Status :: 4 - Beta",
                    "Environment :: Console",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: MIT License",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python",
                    "Topic :: Scientific/Engineering"]
}

if __name__ == '__main__':
    setup(**setup_args)
