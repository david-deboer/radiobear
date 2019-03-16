"""Version definition for radioBEAR."""
from os.path import join as pjoin
import glob

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 2
_version_minor = 4
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
# _version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "Radio-BEAR:  Radio - BErkeley Atmospheric Radiative-transfer"
# Long description will go up on the pypi page
long_description = """
Radio-BEAR
========
Radio - BErkeley Atmospheric Radiative-transfer
To get started using these components in your own software, please go to the
repository README
License
=======
``radio-bear`` is licensed under the terms of the BSD license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.
All trademarks referenced herein are property of their respective holders.
Copyright (c) 2018--, David DeBoer.
"""

NAME = "radio-bear"
MAINTAINER = "David DeBoer"
MAINTAINER_EMAIL = "ddeboer@berkeley.edu"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "https://github.com/david-deboer/radio-bear"
DOWNLOAD_URL = ""
LICENSE = "BSD"
AUTHOR = "David DeBoer"
AUTHOR_EMAIL = "ddeboer@berkeley.edu"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
