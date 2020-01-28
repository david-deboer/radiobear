#! /usr/bin/env python
# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import absolute_import, division, print_function

import glob
import io
from setuptools import setup

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
                 'radiobear.constituents',
                 'radiobear.constituents.clouds', 'radiobear.constituents.co',
                 'radiobear.constituents.h2', 'radiobear.constituents.h2o',
                 'radiobear.constituents.h2s', 'radiobear.constituents.nh3',
                 'radiobear.constituents.ph3'],
    'scripts': glob.glob('scripts/*'),
    'include_package_data': True,
    'install_requires': ["numpy", "matplotlib", "scipy"],
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

# from radiobear import this_radiobear_update  # noqa
