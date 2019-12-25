# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function

# ##local imports
from . import atm_base


class Atmosphere(atm_base.AtmosphereBase):
    def __init__(self, planet, mode='normal', config='config.par', log=None, **kwargs):
        """reads/computes atmospheres.  This computes:
               self.gas
               self.cloud
               self.property
            on the appropriate grid."""
        super(Atmosphere, self).__init__(planet=planet, config=config, log=log)
        self.verbose = False
        if 'verbose' in kwargs.keys():
            self.verbose = kwargs['verbose']
        if self.verbose:
            print('\n---Atmosphere of {}---'.format(planet))

        if self.verbose == 'loud':
            print('Planet ' + self.planet)
            self.config.display()
        self.log.add('\tReading from: ' + self.config.filename, self.verbose)
        self.log.add('\tAtmosphere file:  ' + self.config.gasFile, self.verbose)
        self.log.add('\tCloud file:  ' + self.config.cloudFile, self.verbose)

        self.readGas()
        self.nAtm = len(self.gas[0])
        self.readCloud()
        self.computeProp()
