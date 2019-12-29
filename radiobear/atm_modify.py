# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import sys
import os
import six


def tweakAtm(self):
    """Tweaks the atmosphere data..."""
    if self.batch_mode:
        return 0

    # Import tweakmodule
    __import__(self.config.tweakmodule)
    tweakModule = sys.modules[self.config.tweakmodule]

    # Run module then log
    self.tweakComment, self.gas, self.cloud = tweakModule.modify(self.gas, self.cloud,
                                                                 self.config.C, self.config.Cl)
    if self.verbose:
        print('---tweakComment')
        print(self.tweakComment)
        print('---')
    self.log.add(self.tweakComment, False)
    tf = os.path.join(self.config.path, self.config.tweakmodule + '.py')
    tp = open(tf, 'r')
    dt = tp.read()
    self.log.add('======================' + tf + '=====================', False)
    self.log.add(dt, False)
    self.log.add('====================================================================', False)
    tp.close()


def scaleAtm(self, scale_info='Scratch/scale.dat'):
    """
    This is a built-in tweak module.
    """
    if isinstance(scale_info, six.string_types):
        import alpha
        col, scale_info = alpha.read_scalefile(scale_info)
    else:
        col = scale_info.keys()

    if len(scale_info[col[0]]) != self.nAtm:
        print("Warning - scale file doesn't match atmosphere.  Not applying.")
        return None

    for i in range(self.nAtm):
        for gas in col:
            if gas.lower() != 'p':
                self.gas[self.config.C[gas.upper()]][i] *= scale_info[gas][i]
