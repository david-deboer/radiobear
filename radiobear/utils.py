# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

"""Utilities."""
import os
import numpy as np

commentChars = ['!', '#', '$', '%', '&', '*']
affirmative = [1, '1', 'yes', 'true', 'y', 't']
negative = [0, '0', 'none', 'false', 'n', 'f']

Units = {'Frequency': 1, 'Hz': 1.0, 'kHz': 1.0E3, 'MHz': 1.0E6, 'GHz': 1.0E9,
         'Length': 1, 'm': 1.0, 'km': 1.0E3, 'cm': 1.0E-2, 'AU': 149597870691.0,
         'Pressure': 1, 'bars': 1.0, 'atm': 1.01325,
         'Time': 1, 'sec': 1.0, 'min': 60.0, 'hr': 3600.0, 'day': 86400.0, 'year': 31536000.0,
         'Acceleration': 1, 'mpersec2': 1.0, 'cmpersec2': 0.01}
processingUnits = {'GHz': ['GHz', 'Hz', 'kHz', 'MHz'],
                   'km': ['m', 'cm', 'AU', 'km'],
                   'bars': ['bars', 'atm'],
                   'sec': ['sec', 'min', 'hr', 'day', 'year'],
                   'mpersec2': ['mpersec2', 'cmpersec2']}


def rb_path(add_path=None):
    """Radiobear path."""
    rb_path = os.path.dirname(__file__)
    if add_path is not None:
        rb_path = os.path.join(rb_path, add_path)
    return rb_path


def isanynum(x):
    """Determine if any number type."""
    if isinstance(x, bool) or isinstance(x, dict) or isinstance(x, list):
        return False
    try:
        y = float(x)  # noqa
    except ValueError:
        return False
    return True


def get_location_for_radiobear_setup():
    """Get location for radiobear."""
    return os.path.dirname(__file__)


def proc_unit(supplied_unit):
    """Process a unit."""
    proc_unit = None
    for proc, units in processingUnits.items():
        if supplied_unit in units:
            proc_unit = proc
            break
    return proc_unit


def convert_unit(v, supplied_unit):
    """Convert units."""
    if supplied_unit is None:
        return v
    converted = v
    for proc, units in processingUnits.items():
        if supplied_unit in units:
            converted = v * Units[supplied_unit] / Units[proc]
            break
    return converted


alphaUnit = 'invcm'
atmLayerUnit = 'km'
processingFreqUnit = proc_unit('Hz')
speedOfLight = 2.9979E8     	# m/s
kB = 1.3806503E-23          	# m2 kg s-2 K-1  Boltzmann's constant
hP = 6.626068E-34       	# m2 kg / s	 Planck's constant
T_cmb = 2.725

rfBands = {'HF': [0.003, 0.03], 'VHF': [0.03, 0.3], 'UHF': [0.3, 1.0], 'L': [1.0, 2.0],
           'S': [2.0, 4.0], 'C': [4.0, 8.0], 'X': [8.0, 12.0], 'Ku': [12.0, 18.0],
           'K': [18.0, 26.5], 'Ka': [26.5, 40.0], 'Q': [40.0, 50.0], 'V': [50.0, 75.0],
           'W': [75.0, 110.0]}


def timer(dt):
    """Timer."""
    return dt.seconds + dt.microseconds / 1e6


def r2d(a):
    """Radians to degrees."""
    return a * 180.0 / np.pi


def d2r(a):
    """Degress to radians."""
    return a * np.pi / 180.0


def r2asec(a):
    """Radians to arcseconds."""
    return 3600.0 * r2d(a)


def getRFband(freq, unit='GHz'):
    """Get RF band values."""
    freq = freq * Units[unit] / Units['GHz']
    for bnd in rfBands:
        if rfBands[bnd][0] <= freq < rfBands[bnd][1]:
            return bnd
    return None


def invertDictionary(dic, reverse=False):
    """Invert a dictionary."""
    e = {}
    for d in dic.keys():
        e[dic[d]] = d
    sk = sorted(list(e.keys()), reverse=reverse)
    return e, sk


def ls(directory='Output', tag='dat', show=True, returnList=False):
    """Generate file list for plotTB and writeWavel."""
    filelist = os.listdir(directory)
    files = []
    i = 0
    for fff in filelist:
        show_tag = tag is None or (isinstance(tag, str) and tag in fff)
        if fff[0] != '.' and show_tag:
            files.append(os.path.join(directory, fff))
            if show:
                print('{}:  {}'.format(i, files[i]))
            i += 1
    if returnList:
        return files


def b_type(b):
    """Determine type of b parameter."""
    if isinstance(b, str):
        return b.lower()
    if isinstance(b[0], str):
        return b[0].lower()
    if len(b) > 20:
        return 'image'
    if len(b) > 9:
        return 'profile'
    return 'spectrum'


def get_data_from(line):
    """Get data from line."""
    if line[0] in commentChars or len(line) < 4:
        return None
    try:
        cval = [float(x) for x in line.split()]
    except ValueError:
        cval = None
    return cval


def get_expected_number_of_entries(fp):
    """
    Get expected number of entries per line.

    Get the number of float entries per line and return the most numerous
    length value.  If more than half aren't the maximum value a ValueError
    if raised.
    """
    enoe = {}
    for line in fp:
        cval = get_data_from(line)
        if cval is not None:
            enoe[len(cval)] = enoe.setdefault(len(cval), 0) + 1
    fp.seek(0)
    invdict, sortedkeys = invertDictionary(enoe, reverse=True)

    if len(sortedkeys) > 1:
        if 2 * sum(sortedkeys[1:]) > sortedkeys[0]:
            raise ValueError("Not enough data lines in file: {}".format(sortedkeys))
    return invdict[sortedkeys[0]]


def bang():
    """Show files."""
    from . import plotting
    files = ls(show=False, returnList=True)
    for f in files:
        plotting.plotTB(f, xaxis='wavel', xlog=True, justFreq=True)


# def writeWavel(fn=None, output_file=None, directory='Output'):
#     filename, Tb, f, wavel, b, xlabel, ylabels = readTB(fn=fn, directory=directory)
#     title = filename.split('/')[-1].split('.')[0]
#
#     # Write file
#     title += '_wavel.dat'
#     if output_file is None:
#         output_file = title
#     print('Writing to ', output_file)
#     fp = open(output_file, 'w')
#     for i in range(len(wavel)):
#         s = '%f\t' % (wavel[i])
#         fp.write(s)
#         for j in range(len(b)):
#             s = '%f\t' % (Tb[j][i])
#             fp.write(s)
#         fp.write('\n')
#     fp.close()
#     return n
#
#
# def concatdat(files, directory='Output'):
#     """Given a list of utils.ls indices, returns TB, freq, wavel"""
#     aTB = []
#     bf = []
#     cwavel = []
#     for fil in files:
#         data = readTB(fn=fil)
#         a = list(data[1][0])
#         b = list(data[2])
#         c = list(data[3])
#         aTB += a
#         bf += b
#         cwavel += c
#     sa = np.argsort(np.array(bf))
#     TB = []
#     f = []
#     wavel = []
#     for i in sa:
#         TB.append(aTB[i])
#         f.append(bf[i])
#         wavel.append(cwavel[i])
#     return TB, f, wavel
