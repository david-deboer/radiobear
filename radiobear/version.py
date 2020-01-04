# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

"""Version definition for radioBEAR."""

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 2
_version_minor = 0
_version_micro = 0
_version_extra = ''

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

VERSION = __version__

version_notes = {'1.0': 'First version of radioBEAR from pyplanet.',
                 '1.1': 'General clean-up from transition.',
                 '1.2': 'Make multiple atm/alpha lists.',
                 '2.0': 'Remove state variables and python2.7 support'
                 }
