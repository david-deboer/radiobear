#! /usr/bin/env python
# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import absolute_import, division, print_function

import os
import glob
import io
import json
from setuptools import setup

from Planet import version

with io.open('README.md', 'r', encoding='utf-8') as readme_file:
    readme = readme_file.read()

setup_args = {
    'name': "radioBEAR",
    'description': "radioBEAR:  radio version BErkeley Atmosphere Radiative transfer",
    'long_description': readme,
    'url': "https://github.com/david-deboer/radioBEAR",
    'license': "BSD",
    'author': "David DeBoer",
    'author_email': "ddeboer@berkeley.edu",
    'version': version.VERSION,
    'packages': ['Planet', 'Planet.Constituents', 'Planet.Jupiter', 'Planet.Saturn', 'Planet.Uranus', 'Planet.Neptune'],
    'scripts': glob.glob('scripts/*'),
    'include_package_data': True,
    'install_requires': ["six", "numpy"],
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
