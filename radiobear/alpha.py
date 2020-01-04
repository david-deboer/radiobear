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
        if not self.read_alpha:
            self.setup_formalisms()
        self.saved_fields = ['ordered_constituents', 'alpha_data', 'freqs', 'P']

    def reset_layers(self):
        self.P = None
        self.freqs = None
        self.layers = None
        self.alpha_data = []

    def state(self):
        state_variables.show_state(self)

    def read_scale_file(self, scale_filename):
        return np.load(scale_filename)

    def write_scale_file(self, scale_filename, scale):
        np.save(scale_filename, scale)

    def write_alpha_data(self, save_type='file'):
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

    def read_alpha_data(self, save_type='file'):
        """
        Reads the saved_fields into self.

        Parameters
        ----------
        save_type : str (other are ignored but don't error)
            'file' writes to self.alphafile
            'memory' writes to memory Namespace
        """
        if isinstance(save_type, six.string_types):
            if save_type.lower() == 'file':
                alphain = np.load(self.alphafile)
                for sf in self.saved_fields:
                    setattr(self, sf, alphain[sf])
            elif save_type.lower() == 'memory':
                for sf in self.saved_fields:
                    setattr(self, sf, getattr(self.memory, sf))

    def setup_formalisms(self):
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

    def total_layer_alpha(self, freqs, absorb, scale):
        totalAbsorption = np.zeros_like(freqs)
        if isinstance(scale, float):
            for i in range(len(freqs)):
                totalAbsorption[i] = absorb[i].sum() * scale
            return totalAbsorption

        for i in range(len(freqs)):
            for j in range(absorb.shape[1]):
                constituent = self.ordered_constituents[j]
                new_value = absorb[i, j]
                if constituent in scale.keys():
                    new_value *= scale[constituent]
                totalAbsorption[i] += new_value
            if 'TOTAL' in scale.keys():
                totalAbsorption[i] *= scale['TOTAL']
        return totalAbsorption

    def get_alpha_from_calc(self, freqs, T, P, gas, gas_dict, cloud, cloud_dict, units='invcm'):
        """This gets the total absoprtion coefficient from gas.  It assumes the correct
           frequency units, but maybe should correct that."""
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
        if self.save_alpha:
            self.alpha_data.append(absorb)
        return absorb

    def get_alpha(self, freqs, layer, atm, scale, units='invcm'):
        """This is a wrapper to get the absorption coefficient, either from
           calculating from formalisms or reading from saved data"""

        if self.read_alpha:
            absorb = self.alpha_data[layer]
        else:
            P = atm.gas[atm.config.C['P']][layer]
            T = atm.gas[atm.config.C['T']][layer]
            gas = atm.gas[:, layer]
            cloud = atm.cloud[:, layer]
            absorb = self.get_alpha_from_calc(freqs, T, P, gas, atm.config.C,
                                              cloud, atm.config.Cl, units)
        return self.total_layer_alpha(freqs, absorb, scale)

    def get_layers(self, freqs, atm, scale=False, read_alpha=False, save_alpha=False):
        """
        Compute or read absorption for all layers in atm.

        Note that read_alpha and save_alpha are passed here and overwrite the
        read_alpha/save_alpha state_variables.

        Parameters
        ----------
        freqs : list
            Frequencies to use (in GHz)
        atm : Class Atmosphere
            Atmosphere to use
        scale : dict or float (everything scales to 1.0)
            If dict, needs to be keyed on constituent (or total), if float scales all.
        read_alpha
        save_alpha
        """
        self.freqs = freqs
        self.atm = atm
        self.P = atm.gas[atm.config.C['P']]
        self.read_alpha = read_alpha
        self.save_alpha = save_alpha
        if self.freqs is None and self.save_alpha:
            self.freqs = freqs

        if self.read_alpha:
            self.read_alpha_data(self.read_alpha)
        numLayers = len(atm.gas[0])
        au = utils.alphaUnit
        self.log.add('{} layers'.format(numLayers), self.verbose)
        if isinstance(scale, dict):
            use_varying_scale = True
            for k, v in scale.items():
                if len(v) != numLayers:
                    raise ValueError("Incorrect scale:  N {} vs {}".format(len(v), numLayers))
                if k.lower().startswith('tot'):
                    scale['TOTAL'] = v
        elif isinstance(scale, float):
            use_varying_scale = False
            layer_scale = scale
        else:
            use_varying_scale = False
            layer_scale = 1.0
        layer_alpha = []
        for layer in range(numLayers):
            if use_varying_scale:
                layer_scale = {}
                for k, v in scale.items():
                    layer_scale[k] = v[layer]
            layer_alpha.append(self.get_alpha(freqs, layer, atm, layer_scale, units=au))
        self.layers = np.array(layer_alpha).transpose()
        self.write_alpha_data(save_alpha)
