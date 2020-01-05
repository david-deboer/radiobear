# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import print_function, absolute_import, division

from argparse import Namespace
import numpy as np

from . import utils


def set_b(b, block=[1, 1], **kwargs):
    """
    Process b request.  Sets data_type ('image', 'spectrum', 'profile'), and imSize

       Parameters
       ----------
       b:   list of  pairs -> returns that list
            one pair -> returns that as a list
            float -> returns a full image grid
            'disc' (or 'disk') -> returns [[0.0, 0.0]]
            'stamp:bres:xmin,xmax,ymin,ymax' -> returns grid of postage stamp
            'start:stop:step[<angle] -> string defining line
            'n1,n2,n3[<angle]' -> csv list of b magnitudes
       block:  image block as pair, e.g. [4, 10] is "block 4 of 10"
       kwargs : other args as necessary

       Returns
       -------
       Namespace
           contains b, block, data_type, imSize
       """
    # Deal with strings
    return_value = Namespace(b=None, data_type=None, block=block, imSize=None)
    if isinstance(b, str):
        b = b.lower()
        if b.startswith('dis'):
            return_value.b = [b]
            return_value.data_type = 'spectrum'
        if b.startswith('stamp'):
            bres = float(b.split(':')[1])
            bext = [float(x) for x in b.split(':')[2].split(',')]
            return_value.b = []
            for x in np.arange(bext[0], bext[1] + bres / 2.0, bres):
                for y in np.arange(bext[2], bext[3] + bres / 2.0, bres):
                    return_value.b.append([y, x])
            xbr = len(np.arange(bext[2], bext[3] + bres / 2.0, bres))
            return_value.imSize = [xbr, len(b) / xbr]
            return_value.data_type = 'image'
        else:
            b = b.split('<')
            angle_b = 0.0 if len(b) == 1 else utils.d2r(float(b[1]))
            if ',' in b[0]:
                mag_b = [float(x) for x in b[0].split(',')]
            elif ':' in b[0]:
                mag_b = [float(x) for x in b[0].split(':')]
                mag_b = np.arange(mag_b[0], mag_b[1] + mag_b[2] / 2.0, mag_b[2])
            ab = kwargs['Rpol'] / kwargs['Req']
            rab = ab / np.sqrt(np.power(np.sin(angle_b), 2.0) + np.power(ab * np.cos(angle_b), 2.0))
            return_value.b = []
            for v in mag_b:
                if v < 0.99 * rab:
                    return_value.b.append([v * np.cos(angle_b), v * np.sin(angle_b)])
            return_value.data_type = 'profile'
        return return_value
    if isinstance(b, float):  # this generates a grid at that spacing and blocking
        grid = -1.0 * np.flipud(np.arange(b, 1.5 + b, b))
        grid = np.concatenate((grid, np.arange(0.0, 1.5 + b, b)))
        # get blocks
        bsplit = len(grid) / abs(block[1])
        lastRow = block[0] / abs(block[1])
        if abs(block[1]) == 1:
            lastRow = 0
        return_value.b = []
        for i in range(int(bsplit + lastRow)):
            ii = i + int((block[0] - 1) * bsplit)
            vrow = grid[ii]
            for vcol in grid:
                return_value.b.append([vcol, vrow])
        return_value.imSize = [len(grid), len(b) / len(grid)]
        return_value.data_type = 'image'
        return return_value
    shape_b = np.shape(b)
    if len(shape_b) == 1:
        return_value.data_type = 'spectrum'
        return_value.b = [b]
    else:
        return_value.data_type = 'spectrum' if shape_b[0] < 5 else 'profile'
        return_value.b = b
    return return_value


def set_freq(freqs, freqUnit='GHz'):
    """
    Process frequency request.
    Return a list converted from freqUnit to processingFreqUnit and reassigns freqUnit procUnit.

    Parameters
    ----------
    freqs:  list -> returns that list
            csv list of values -> converts to list
            float/int -> returns that one value as a list
            '<start>:<stop>:<step>' -> returns arange of that
            '<start>;<stop>;<nstep>' -> returns logspace of that
            '<filename>' -> returns loadtxt of that file
    freqUnit : str
        Frequency unit of supplied freqs

    Returns
    -------
    list
        Frequency list
    str
        Frequency unit
    """
    if isinstance(freqs, list):
        pass
    elif isinstance(freqs, np.ndarray):
        freqs = list(freqs)
    elif isinstance(freqs, str) and ',' in freqs:
        freqs = [float(x) for x in freqs.split(',')]
    elif isinstance(freqs, (float, int)):
        freqs = [float(freqs)]
    elif isinstance(freqs, str) and ':' in freqs:
        fstart, fstop, fstep = [float(x) for x in freqs.split(':')]
        freqs = list(np.arange(fstart, fstop + fstep / 2.0, fstep))
    elif isinstance(freqs, str) and ';' in freqs:
        fstart, fstop, nstep = [float(x) for x in freqs.split(';')]
        freqs = list(np.logspace(np.log10(fstart), np.log10(fstop), nstep))
    elif isinstance(freqs, str):
        freqs = list(np.loadtxt(freqs))
    else:
        raise ValueError('Invalid format for frequency request')

    for i in range(len(freqs)):
        freqs[i] = utils.convert_unit(freqs[i], freqUnit)

    return freqs, freqUnit
