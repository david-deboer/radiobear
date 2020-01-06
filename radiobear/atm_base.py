# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
import numpy as np
import os
import sys

# ##local imports
from . import utils
from . import chemistry
from . import logging


class AtmosphereBase:
    def __init__(self, planet, config=None, log=None, **kwargs):
        """
        Atmosphere base class
        """
        self.planet = planet.capitalize()
        self.log = logging.setup(log)

        if config is None or isinstance(config, str):
            from . import config as pcfg
            config = os.path.join(self.planet, config)
            config = pcfg.planetConfig(self.planet, configFile=config)
            config.update_config(**kwargs)
        self.config = config
        self.chem = None
        self.gas = []
        self.nGas = 0
        self.nConstituent = 0
        self.cloud = []
        self.nCloud = 0
        self.property = []
        if config.path not in sys.path:
            sys.path.insert(0, config.path)

    def readGas(self):
        """Reads gas profile file as self.gas"""
        gasFile = os.path.join(self.config.path, self.config.gasFile[self.idnum])

        self.gas = []
        Cdict = self.config.C
        value_sorted_Cdict_keys = sorted(Cdict, key=lambda k: Cdict[k])
        for k in value_sorted_Cdict_keys:
            self.gas.append([])
        self.nConstituent = len(Cdict)
        with open(gasFile, "r") as fp:
            expected_number = utils.get_expected_number_of_entries(fp)
            for line in fp:
                cval = utils.get_data_from(line)
                if cval is None or len(cval) != expected_number:
                    continue
                for n in range(self.nConstituent):
                    if n < len(cval):
                        self.gas[n].append(cval[n])
                    else:
                        self.gas[n].append(0.0)
        self.nGas = len(self.gas[0])
        # ##Redo the constituent dictionary for the self.gas index positions
        nid, sk = utils.invertDictionary(Cdict)
        for i, k in enumerate(sk):
            Cdict[nid[k]] = i
        self.config.C = Cdict

        self.gas = np.array(self.gas)
        # ## Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.gas[self.config.C['P']]) > 0.0)
        if not monotonic:
            self.gas = np.fliplr(self.gas)
            monotonic = np.all(np.diff(self.gas[self.config.C['P']]) > 0.0)
        if not monotonic:
            raise ValueError("Pressure not monotonically increasing in {}.".format(gasFile))

        # ## Renormalize so that deepest z is 0 and set DZ
        self._renorm_z('gas')

    def readCloud(self):
        """Reads in cloud data if we have it..."""

        Cldict = self.config.Cl
        cloudFile = os.path.join(self.config.path, self.config.cloudFile[self.idnum])

        self.cloud = []
        value_sorted_Cldict_keys = sorted(Cldict, key=lambda k: Cldict[k])
        for k in value_sorted_Cldict_keys:
            self.cloud.append([])
        self.nParticulate = len(Cldict.keys())
        with open(cloudFile, "r") as fp:
            expected_number = utils.get_expected_number_of_entries(fp)
            for line in fp:
                cval = utils.get_data_from(line)
                if cval is None or len(cval) != expected_number:
                    continue
                for n in range(self.nParticulate):  # Initialize all of the particulates to 0.0
                    if n < len(cval):
                        self.cloud[n].append(cval[n])
                    else:
                        self.cloud[n].append(0.0)
        self.nCloud = len(self.cloud[0])
        # ##Redo the particulate dictionary for the self.cloud index positions
        nid, sk = utils.invertDictionary(Cldict)
        for i, k in enumerate(sk):
            Cldict[nid[k]] = i
        self.config.Cl = Cldict

        self.cloud = np.array(self.cloud)
        # ## Check that P is monotonically increasing
        monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']]) > 0.0)
        if not monotonic:
            self.cloud = np.fliplr(self.cloud)
            monotonic = np.all(np.diff(self.cloud[self.config.Cl['P']]) > 0.0)
        if not monotonic:
            raise ValueError("Pressure not monotonically increasing in {}".format(cloudFile))

        # ## Renormalize so that deepest z is 0 and set DZ
        self._renorm_z('cloud')

    def computeProp(self):
        """This module computes derived atmospheric properties (makes self.property)"""
        self.chem = {}
        for key in self.config.C:
            if key in ['P', 'T', 'Z', 'DZ']:
                continue
            self.chem[key] = chemistry.ConstituentProperties(key)

        # nAtm = len(self.gas[self.config.C['P']])
        self.property = []
        for op in self.config.LP:
            self.property.append([])
        zOffset = 0.0
        iOffset = 0
        psep = 1.0E6
        for i, zv in enumerate(self.gas[self.config.C['Z']]):  # find the nearest z value at p_ref
            P = self.gas[self.config.C['P']][i]
            if abs(P - self.config.p_ref) < psep:
                psep = abs(P - self.config.p_ref)
                iOffset = i
        zOffset = self.gas[self.config.C['Z']][iOffset]
        z_at_p_ref = self.config.Req

        for i, zv in enumerate(self.gas[self.config.C['Z']]):
            T = self.gas[self.config.C['T']][i]
            P = self.gas[self.config.C['P']][i]
            self.property[self.config.LP['P']].append(P)
            self.property[self.config.LP['Z']].append(zv)
            rr = z_at_p_ref + zv - zOffset
            # note that this is the "actual"z along equator referenced to planet center (aka radius)
            self.property[self.config.LP['R']].append(rr)
            # ##set mean amu
            amulyr = 0.0
            for key in self.chem:
                amulyr += self.chem[key].amu * self.gas[self.config.C[key]][i]
            self.property[self.config.LP['AMU']].append(amulyr)
            # ##set GM pre-calc (normalized further down) and get lapse rate
            if not i:
                self.property[self.config.LP['GM']].append(0.0)
                self.property[self.config.LP['LAPSE']].append(0.0)
                self.property[self.config.LP['LAPSEP']].append(0.0)
            else:
                rho = (amulyr * P) / (chemistry.R * T)
                dr = abs(zv - self.gas[self.config.C['Z']][i - 1])
                dV = 4.0 * np.pi * (rr**2) * dr
                dM = 1.0e11 * rho * dV
                GdM = self.property[self.config.LP['GM']][i - 1] + chemistry.GravConst * dM
                # in km3/s2
                # mass added as you make way into atmosphere by radius r (times G)
                self.property[self.config.LP['GM']].append(GdM)
                dT = abs(T - self.gas[self.config.C['T']][i - 1])
                dP = abs(P - self.gas[self.config.C['P']][i - 1])
                self.property[self.config.LP['LAPSE']].append(dT / dr)
                self.property[self.config.LP['LAPSEP']].append(dT / dP)
            # ##set refractivity and index of refraction
            refrlyr = 0.0
            for key in self.chem:
                refrlyr += self.chem[key].refractivity(T=T) * self.gas[self.config.C[key]][i]
            refrlyr = refrlyr * P * (293.0 / T)
            self.property[self.config.LP['REFR']].append(refrlyr)
            nlyr = refrlyr / 1.0E6 + 1.0
            self.property[self.config.LP['N']].append(nlyr)

        # ##Now need to normalize GM to planet and calculate scale height (H)
        GMnorm = self.property[self.config.LP['GM']][iOffset]  # G*(Mass added by p_ref)
        for i, mv in enumerate(self.property[self.config.LP['GM']]):
            gm = self.config.GM_ref - (mv - GMnorm)
            self.property[self.config.LP['GM']][i] = gm
            little_g = gm / self.property[self.config.LP['R']][i]**2
            m_bar = self.property[self.config.LP['AMU']][i]
            T = self.gas[self.config.C['T']][i]
            self.property[self.config.LP['H']].append((chemistry.R * T) /
                                                      (little_g * m_bar) / 1000.0)
            self.property[self.config.LP['g']].append(little_g)
        self.property = np.array(self.property)

    def is_present(self, c, tiny=1.0E-30):
        """This checks to see if a constituent is there and sets 0.0 or negative values to tiny.
           This is generally for log plotting."""
        v = [tiny if x <= tiny else x for x in c]
        present = bool(len(np.where(np.array(v) > tiny)[0]))
        return present, v

    def _renorm_z(self, gctype):
        if gctype == 'gas':
            ind = self.config.C
            atm_gc = self.gas
        else:
            ind = self.config.Cl
            atm_gc = self.cloud
        zDeep = atm_gc[ind['Z']][-1]
        for i in range(len(atm_gc[ind['Z']])):
            atm_gc[ind['Z']][i] -= zDeep

        # put in DZ
        dz = np.abs(np.diff(atm_gc[ind['Z']])) * 1.0E5  # convert from km to cm - no unit checking!!
        atm_gc[ind['DZ']] = np.append(np.array([0.0]), dz)
