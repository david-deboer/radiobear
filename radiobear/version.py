# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

"""Version definition for radioBEAR."""

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 1
_version_minor = 1
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
