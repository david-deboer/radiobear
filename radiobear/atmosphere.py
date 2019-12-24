# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import numpy as np
import sys
import os
import six

# ##local imports
from . import utils
from . import regrid
from . import state_variables
from . import atm_base


class Atmosphere(atm_base.AtmosphereBase):
    def __init__(self, planet, mode='normal', config='config.par', log=None, **kwargs):
        """reads/computes atmospheres.  This computes:
               self.gas
               self.cloud
               self.property
            on the appropriate grid."""
        super(Atmosphere, self).__init__(planet=planet, config=config, log=log)
        kwargs = state_variables.init_state_variables(mode, **kwargs)
        self.state_vars = kwargs.keys()
        state_variables.set_state(self, set_mode='init', **kwargs)
        if self.verbose:
            print('\n---Atmosphere of {}---'.format(planet))

        # ##Create function dictionaries
        self.gasGen = {}
        self.gasGen['read'] = self.readGas
        self.cloudGen = {}
        self.cloudGen['read'] = self.readCloud
        self.propGen = {}
        self.propGen['compute'] = self.computeProp

        if self.verbose == 'loud':
            print('Planet ' + self.planet)
            self.config.display()
        if self.config.gasType == 'read':  # this assumes that cloudType is then also 'read'
            self.log.add('\tReading from: ' + self.config.filename, self.verbose)
            self.log.add('\tAtmosphere file:  ' + self.config.gasFile, self.verbose)
            self.log.add('\tCloud file:  ' + self.config.cloudFile, self.verbose)

    def state(self):
        state_variables.show_state(self)

    def run(self, Pmin=None, Pmax=None, regridType=None, gasType=None, cloudType=None,
            otherType=None, tweak=True):
        """This is the standard pipeline"""
        # ##Set run defaults
        if Pmin is None:
            Pmin = self.config.pmin
        if Pmax is None:
            Pmax = self.config.pmax
        if regridType is None:
            regridType = self.config.regridType
        if gasType is None:
            gasType = self.config.gasType
        if cloudType is None:
            cloudType = self.config.cloudType
        if otherType is None:
            otherType = self.config.otherType

        # ## Generate gas profile
        if gasType not in self.gasGen.keys():
            print('Error:  No such gasType: ', gasType)
            return 0
        else:
            self.gasGen[gasType]()
        self.nAtm = len(self.gas[0])

        if self.batch_mode:
            return self.nAtm

        # ## Generate cloud profile
        if cloudType not in self.cloudGen.keys():
            print('Error:  No such cloudType: ', cloudType)
            return 0
        else:
            self.cloudGen[cloudType]()

        # ## Put onto common grid
        if self.verbose:
            print("Regrid:  {}".format(regridType))
        regridded = regrid.regrid(self, regridType=regridType, Pmin=Pmin, Pmax=Pmax)  # noqa
        self.nAtm = len(self.gas[0])

        if tweak:  # This loads/calls the module as given in the config.par tweakmodule parameter
            self.tweakAtm()

        # ## Compute other parameters that are needed
        if self.config.otherType not in self.propGen.keys():
            print('Error:  no such otherTpe: ', otherType)
            return 0
        else:
            self.propGen[otherType]()

        angularDiameter = 2.0 * np.arctan(self.property[self.config.LP['R']][0] /
                                          self.config.distance)
        if self.verbose == 'loud':
            print('angular radius = {} arcsec'.format(utils.r2asec(angularDiameter / 2.0)))

        return self.nAtm

    def tweakAtm(self):
        """Tweaks the atmosphere data..."""
        if self.batch_mode:
            return 0

        # Import tweakmodule
        sys.path.append(self.config.path)
        try:
            __import__(self.config.tweakmodule)
            tweakModule = sys.modules[self.config.tweakmodule]
        except SyntaxError:
            self.log.add("Syntax Error:  check " + self.config.tweakmodule, True)
            raise ValueError("Error in tweakmodule")

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
