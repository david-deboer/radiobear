# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

import numpy as np
import six
import copy
from radiobear import utils


class Data:
    """
    This holds the data that one may wish to use as output.  Data are stored as numpy arrays
    """

    allowed_parameters = ['f', 'freqUnit', 'b', 'Tb', 'header', 'start', 'stop',
                          'log', 'type', 'logfile']

    def __init__(self):
        for a in self.allowed_parameters:
            setattr(self, a, None)

    def __repr__(self):
        s = ''
        for i, b in enumerate(self.b):
            bstr = 'b = {}:'.format(b)
            s += bstr + '\n f = '
            f = ['{:6.1f}'.format(x) for x in self.f]
            s += ' '.join(f) + '  GHz\nTb = '
            T = ['{:6.1f}'.format(x) for x in self.Tb[i]]
            s += ' '.join(T) + '  K\n'
        return s

    def set(self, par, val):
        if par not in self.allowed_parameters:
            print("{} not in valid data return list.".format(par))
            return
        val = copy.copy(val)
        if par == 'b' and utils.b_type(val).startswith('dis'):
            self.b = ['disc']
        elif isinstance(val, list):
            setattr(self, par, np.asarray(val, dtype=np.float32))
        else:
            setattr(self, par, val)

    def show(self, include=['header', 'start', 'stop', 'f', 'b', 'Tb', 'log'], indent=1):
        tab = indent * '\t'
        for v in include:
            if v == 'log':
                self.show_log(indent=indent)
                continue
            if v not in self.allowed_parameters:
                continue
            if v == 'header':
                self.show_header(indent=indent)
            elif v == 'f':
                print("{}freq: {} {}".format(tab, self.f, self.freqUnit))
            else:
                print("{}{}:  {}".format(tab, v, getattr(self, v)))

    def show_header(self, indent=1):
        tab = indent * '\t'
        print('{}Header'.format(tab))
        tab1 = (indent + 1) * '\t'
        for k, v in six.iteritems(self.header):
            print("{}{:20s}     {}".format(tab1, k, v))

    def show_log(self, indent=1):
        tab = indent * '\t'
        if self.log is not None:
            self.log.close()
        print('{}Log:  {}'.format(tab, self.logfile))
        with open(self.logfile) as fp:
            for line in fp:
                print(line.strip())
        print("-----------------------------------")
