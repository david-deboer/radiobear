# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function

# ##local imports
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
        state_variables.init_state_variables(self, mode, **kwargs)
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

    def simple(self):
        self.readGas()
        self.nAtm = len(self.gas[0])
        self.readCloud()
        self.computeProp()

    def std(self, Pmin=None, Pmax=None, regridType=None, gasType=None, cloudType=None,
            otherType=None, tweak=True):
        from . import regrid
        from . import atm_modify
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
            atm_modify.tweakAtm(self)

        # ## Compute other parameters that are needed
        if self.config.otherType not in self.propGen.keys():
            print('Error:  no such otherTpe: ', otherType)
            return 0
        else:
            self.propGen[otherType]()

        return self.nAtm
