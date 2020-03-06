# -*- mode: python; coding: utf-8 -*-
"""Atmosphere class on top of base."""
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

from . import atm_base


class Atmosphere(atm_base.AtmosphereBase):
    """Atmosphere class on top of base."""

    def __init__(self, planet, idnum=0, config='config.par', log=None, verbose=False, **kwargs):
        """
        Reads/computes the atmosphere to be used.

        This computes:
               self.gas
               self.cloud
               self.property
            on the appropriate grid and for the supplied idnum.

        Parameters
        ----------
        planet : str
            Planet name
        idnum : int
            For read, index number of the gasFile/cloudFile to be used.
        config : str or class
            Configuration to use
        log : None or str
            Log setup
        """
        super().__init__(planet=planet, config=config, log=log, **kwargs)
        self.verbose = verbose
        if self.verbose:
            print('\n---Atmosphere of {}---'.format(planet))

        # ##Create function dictionaries
        self.gasGen = {}
        self.gasGen['read'] = self.readGas
        self.cloudGen = {}
        self.cloudGen['read'] = self.readCloud
        self.propGen = {}
        self.propGen['compute'] = self.computeProp

        self.idnum = idnum
        if self.verbose == 'loud':
            print('Planet ' + self.planet)
            print(self.config)
        if self.config.gasType == 'read':  # this assumes that cloudType is then also 'read'
            self.log.add('\tReading from: ' + self.config.filename, self.verbose)
            self.log.add('\tAtmosphere file:  ' + str(self.config.gasFile[idnum]), self.verbose)
            self.log.add('\tCloud file:  ' + str(self.config.cloudFile[idnum]), self.verbose)

    def simple(self, **kwwargs):
        """
        Compute simple pipeline, that uses the files as-is.

        kwargs just catches any values that may needed for std below.
        """
        self.readGas()
        self.readCloud()
        self.computeProp()
        self.nAtm = len(self.gas[0])

    def std(self, Pmin=None, Pmax=None, regridType=None, gasType=None, cloudType=None,
            otherType=None, tweak=True):
        """Compute the standard pipeline."""
        from . import regrid
        from . import atm_modify
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
