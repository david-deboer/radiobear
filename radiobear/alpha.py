# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import os
import sys
import six
import numpy as np
from . import utils
from . import config as pcfg
from . import state_variables
from . import logging


class Alpha:
    def __init__(self, idnum=0, mode='normal', config=None, log=None, **kwargs):
        """
        Reads in absorption formalisms and computes layer absorption.  Note that they are all in GHz

        Parameters
        ----------
        planet : str
            Planet name
        idnum : int
            For read, index number of the gasFile/cloudFile to be used.
        mode : str
            Mode of atmosphere use.  Look under "state_variables."
        config : str or class
            Configuration to use
        log : None or str
            Log setup
        **kwargs
        """

        state_variables.init_state_variables(self, mode, **kwargs)
        self.log = logging.setup(log)
        self.constituentsAreAt = os.path.join(os.path.dirname(__file__), 'Constituents')
        self.idnum = idnum
        self.reset_layers()

        # get config
        if config is None or isinstance(config, str):
            config = pcfg.planetConfig(self.planet, configFile=config, log=self.log)
        self.config = config

        self.alphafile = os.path.join(self.scratch_directory, 'alpha{:04d}.npz'.format(self.idnum))
        if self.read_alpha:
            self.existing_alpha_setup()
        else:
            self.formalisms()
        self.saved_fields = ['ordered_constituents', 'alpha_data', 'freqs', 'P']

        # copy config back into other_dict as needed
        other_to_copy = {}
        other_to_copy['h2'] = ['h2state', 'h2newset']
        other_to_copy['clouds'] = ['water_p', 'ice_p', 'nh4sh_p', 'nh3ice_p', 'h2sice_p', 'ch4_p']
        self.other_dict = {}
        for absorber in self.ordered_constituents:
            self.other_dict[absorber] = {}
            if absorber in other_to_copy.keys():
                for oc in other_to_copy[absorber]:
                    self.other_dict[absorber][oc] = getattr(self.config, oc)

    def reset_layers(self):
        self.P = None
        self.freqs = None
        self.layers = None
        self.alpha_data = []

    def get_layers(self, freqs, atm, read_alpha=False, save_alpha=False):
        self.freqs = freqs
        self.atm = atm
        self.P = atm.gas[atm.config.C['P']]
        numLayers = len(atm.gas[0])
        layerAlp = []
        self.log.add('{} layers'.format(numLayers), self.verbose)
        for layer in range(numLayers):
            layerAlp.append(self.get_alpha(freqs, layer, atm, units=utils.alphaUnit))
            if self.verbose == 'loud':
                print('\r\tAbsorption in layer {}   '.format(layer + 1), end='')
        self.layers = np.array(layerAlp).transpose()
        self.write_alpha(save_alpha)

    def read_scale(self, scale_type):
        """
        """
        if isinstance(scale_type, six.string_types):
            try:
                self.scale = float(scale_type)
            except ValueError:
                self.scale = np.loadtxt(scale_type)
        elif scale_type:
            self.scale = getattr(scale_type, 'scale')

    def write_alpha(self, save_type='file'):
        """
        Write out the saved_fields.

        Parameters
        ----------
        save_type : str (other are ignored but don't error)
            'file' writes to self.alphafile
            'memory' writes to memory Namespace
        """
        if isinstance(save_type, six.string_types):
            self.alpha_data = np.array(self.alpha_data)
            if save_type.lower() == 'file':
                np.savez(self.alphafile,
                         ordered_constituents=self.ordered_constituents,
                         alpha_data=self.alpha_data,
                         freqs=self.freqs, P=self.P)
            elif save_type.lower() == 'memory':
                for sf in self.saved_fields:
                    setattr(self.memory, sf, getattr(self, sf))

    def read_alpha(self, save_type='file'):
        """
        Reads the saved_fields into self.

        Parameters
        ----------
        save_type : str, bool, None, Namespace, Class
                    if str or True, assumes read from file self.alphafile
                    if Namespace or Class, reads from that
                    if not bool(save_type), does nothing
                    otherwise, raises error
        """
        if isinstance(save_type, bool) and save_type:
            save_type = 'file'
        if isinstance(save_type, str):
            alphain = np.load(self.alphafile)
            for sf in self.saved_fields:
                setattr(self, sf, alphain[sf])
        elif save_type:
            for sf in self.saved_fields:
                setattr(self, sf, getattr(save_type, sf))

    def formalisms(self):
        # Get possible constituents
        s = 'Reading in absorption modules from ' + self.constituentsAreAt + '\n'
        self.log.add(s, self.verbose)
        # Import used ones - note this dynamically imports the absorption modules.
        self.constituent = {}
        self.absorptionModule = {}
        self.truncate_strength = {}
        self.truncate_freq = {}
        for c in self.config.constituent_alpha:
            absorber = self.config.constituent_alpha[c]
            if absorber is None:
                continue
            constituentPath = os.path.join(self.constituentsAreAt, c)
            sys.path.append(constituentPath)
            try:
                __import__(absorber)
                self.absorptionModule[c] = sys.modules[absorber]
                self.constituent[c] = absorber
            except ImportError:
                s = "WARNING:  CAN'T LOAD " + absorber + '.  '
                print(s * 3)
                self.log.add("Can't load " + absorber, True)
            self.truncate_freq[c] = None
            self.truncate_strength[c] = None
            if c in self.config.truncate_method.keys() and\
                    self.config.truncate_method[c] is not None:
                trunc_methods = self.config.truncate_method[c].split(",")
                for tm in trunc_methods:
                    getattr(self, tm)[c] = getattr(self.config, tm)[c]
        self.ordered_constituents = sorted(self.constituent.keys())
        self.log.add('Using modules:', self.verbose)
        for c in self.constituent:
            s = "\t{}: {} \t".format(c, self.constituent[c])
            if self.truncate_strength[c] is not None:
                s += "(truncate_strength {})\t".format(self.truncate_strength[c])
            if self.truncate_freq[c] is not None:
                s += "(truncate_freq {})".format(self.truncate_freq[c])
            self.log.add(s, self.verbose)

    def get_alpha(self, freqs, layer, atm, units='invcm'):
        """This is a wrapper to get the absorption coefficient, either from
           calculating from formalisms or reading from file"""
        if self.freqs is None and self.save_alpha:
            self.freqs = freqs
        if self.read_alpha:
            if len(self.alpha_data) != len(atm.gas[0]):
                raise ValueError("Absorption and atmosphere don't agree")
        if self.read_alpha:
            return self.get_alpha_from_file(freqs, layer, units)
        else:
            P = atm.gas[atm.config.C['P']][layer]
            T = atm.gas[atm.config.C['T']][layer]
            gas = atm.gas[:, layer]
            cloud = atm.cloud[:, layer]
            return self.get_alpha_from_calc(freqs, T, P, gas, atm.config.C,
                                            cloud, atm.config.Cl, units)

    def get_alpha_from_file(self, freqs, layer, units='invcm'):
        totalAbsorption = np.zeros_like(freqs)
        if isinstance(self.scale_by, float):
            for i in range(len(freqs)):
                totalAbsorption[i] = self.alpha_data[layer, i].sum() * self.scale_by
            return totalAbsorption

        for i in range(len(freqs)):
            for j in range(self.alpha_data.shape[2] - 1):
                constituent = self.ordered_constituents[j]
                new_value = self.alpha_data[layer, i, j]
                if constituent in self.scale_constituent_columns:
                    new_value *= self.scale_constituent_values[constituent][layer]
                totalAbsorption[i] += new_value
        return totalAbsorption

    def get_alpha_from_calc(self, freqs, T, P, gas, gas_dict, cloud, cloud_dict, units='invcm'):
        """This gets the total absoprtion coefficient from gas.  It assumes the correct
           frequency units, but maybe should correct that.
           Returns total absorption at that layer."""
        absorb = []
        print_meta = self.verbose == 'loud'
        for c in self.ordered_constituents:
            path = os.path.join(self.constituentsAreAt, c)
            if c.lower().startswith('cloud'):
                X = cloud
                D = cloud_dict
            else:
                X = gas
                D = gas_dict
            absorb.append(self.absorptionModule[c].alpha(freqs, T, P, X, D, self.other_dict[c],
                          truncate_freq=self.truncate_freq[c],
                          truncate_strength=self.truncate_strength[c],
                          units=units, path=path, verbose=print_meta))
        absorb = np.array(absorb)
        absorb = absorb.transpose()
        totalAbsorption = np.zeros_like(freqs)
        for i in range(len(freqs)):
            totalAbsorption[i] = absorb[i].sum()
        if self.save_alpha:
            this_layer = []
            for this_freq in absorb:
                this_layer.append(this_freq)
            self.alpha_data.append(this_layer)
        del absorb
        return totalAbsorption

    def state(self):
        state_variables.show_state(self)
