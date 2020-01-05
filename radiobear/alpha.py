# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import os
import sys
import numpy as np
from . import utils
from . import logging


class Alpha:
    def __init__(self, idnum=0, config=None, log=None, load_formal=True, verbose=True, **kwargs):
        """
        Reads in absorption formalisms and computes layer absorption.  Note that they are all in GHz

        Parameters
        ----------
        idnum : int
            For read, index number of the gasFile/cloudFile to be used.
        config : str or class
            Configuration to use
        log : None or str
            Log setup
        load_formal : bool
            Flag to load in the absorption modules
        """
        self.verbose = verbose
        self.log = logging.setup(log)
        if config is None or isinstance(config, str):
            from . import config as pcfg
            config = pcfg.planetConfig('x', configFile=config)
            config.update_config(**kwargs)
        self.config = config
        self.constituentsAreAt = os.path.join(os.path.dirname(__file__), 'Constituents')
        self.idnum = idnum
        self.reset_layers()
        self.alphafile = os.path.join(self.config.scratch_directory,
                                      'alpha{:04d}.npz'.format(self.idnum))
        if load_formal:
            self.setup_formalisms()
        self.saved_fields = ['ordered_constituents', 'alpha_data', 'freqs', 'P']

    def reset_layers(self):
        self.P = None
        self.freqs = None
        self.layers = None
        self.alpha_data = []

    def save_alpha_data(self, save_type='file'):
        """
        Write out the saved_fields.

        Parameters
        ----------
        save_type : str (other are ignored but don't error)
            'file' writes to self.alphafile
            'memory' writes to memory Namespace
        """
        self.alpha_data = np.array(self.alpha_data)
        if save_type == 'file':
            np.savez(self.alphafile,
                     ordered_constituents=self.ordered_constituents,
                     alpha_data=self.alpha_data,
                     freqs=self.freqs, P=self.P)
        elif save_type == 'memory':
            for sf in self.saved_fields:
                setattr(self.memory, sf, getattr(self, sf))

    def read_alpha_data(self, save_type):
        """
        Reads the saved_fields into self.

        Parameters
        ----------
        save_type : str (other are ignored but don't error)
            'file' writes to self.alphafile
            'memory' writes to memory Namespace
        """
        if save_type == 'file':
            alphain = np.load(self.alphafile)
            for sf in self.saved_fields:
                setattr(self, sf, alphain[sf])
        elif save_type == 'memory':
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

    def total_layer_alpha(self, absorb, lscale):
        """
        Takes the layer absorption profile and scale-sums it.

        Parameters
        ----------
        """
        totalAbsorption = np.zeros_like(self.freqs)
        Nfreq = len(self.freqs)
        if self._save_alpha_memfil:
            save_absorb = np.zeros_like(absorb)
        if isinstance(lscale, float):
            for i in range(Nfreq):
                totalAbsorption[i] = absorb[i].sum() * lscale
                if self._save_alpha_memfil:
                    for j in absorb[i]:
                        save_absorb[i, j] = absorb[i, j] * lscale
            if self._save_alpha_memfil:
                self.alpha_data.append(save_absorb)
            return totalAbsorption

        for i in range(Nfreq):
            for j in range(absorb.shape[1]):
                constituent = self.ordered_constituents[j]
                new_value = absorb[i, j]
                if constituent in lscale.keys():
                    new_value *= lscale[constituent]
                if self._save_alpha_memfil:
                    save_absorb[i, j] = new_value
                totalAbsorption[i] += new_value
        if self._save_alpha_memfil:
            self.alpha_data.append(save_absorb)
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
        return absorb

    def get_single_layer(self, freqs, layer, atm, lscale, units='invcm'):
        """This is a wrapper to get the absorption coefficient, either from
           calculating from formalisms or reading from saved data"""

        if self._get_alpha_memfil:
            absorb = self.alpha_data[layer]
        else:
            P = atm.gas[atm.config.C['P']][layer]
            T = atm.gas[atm.config.C['T']][layer]
            gas = atm.gas[:, layer]
            cloud = atm.cloud[:, layer]
            absorb = self.get_alpha_from_calc(freqs, T, P, gas, atm.config.C,
                                              cloud, atm.config.Cl, units)
        return self.total_layer_alpha(absorb=absorb, lscale=lscale)

    def get_layer_scale(self, scale, N):
        layer_scale = []
        if isinstance(scale, dict):
            nvalid = 0
            for k, v in scale.items():
                if k not in self.ordered_constituents:
                    if self.verbose:
                        print("{} not found as constituent for alpha, so ignoring.".format(k))
                else:
                    nvalid += 1
                if not nvalid:
                    raise ValueError("No valid scale constituents found.")
                if len(v) != N:
                    raise ValueError("Incorrect scale for {}:  N {} vs {}".format(k, len(v), N))
            for lyr in range(N):
                layer_scale.append({})
                for k, v in scale.items():
                    layer_scale[lyr][k] = v[lyr]
        elif isinstance(scale, (list, np.array)):
            if len(scale) != N:
                raise ValueError("Incorrect number of scale layers.")
            layer_scale = scale
        else:
            if isinstance(scale, (float, int)):
                scale_val = float(scale)
            else:
                scale_val = 1.0
            for lyr in range(N):
                layer_scale.append(scale_val)
        return layer_scale

    def get_layers(self, freqs, atm, scale=False, get_alpha='calc', save_alpha='none'):
        """
        Compute or read absorption for all layers in atm.

        Parameters
        ----------
        freqs : list
            Frequencies to use (in GHz)
        atm : Class Atmosphere
            Atmosphere to use
        scale : dict, list or float/int.
            If dict, needs to be keyed on constituent.
            If list, scales total at each layer.
            If float/int, scales total for all layers.
            Other, scales at 1.0
        get_alpha : str
            How to get the absorption:  'file', 'memory', 'calc'
        save_alpha : str
            If/how to save the absoprtion:  'file', 'memory', 'none'
        """
        self.freqs = freqs
        self.atm = atm
        self.P = atm.gas[atm.config.C['P']]
        self.get_alpha = get_alpha
        self._get_alpha_memfil = get_alpha in ['memory', 'file']
        self.save_alpha = save_alpha
        self._save_alpha_memfil = save_alpha in ['memory', 'file']

        self.read_alpha_data(self.get_alpha)
        numLayers = len(atm.gas[0])
        au = utils.alphaUnit
        self.log.add('{} layers'.format(numLayers), self.verbose)
        lscale = self.get_layer_scale(scale, numLayers)
        layer_alpha = []
        for layer in range(numLayers):
            layer_alpha.append(self.get_single_layer(freqs, layer, atm, lscale[layer], units=au))
        self.layers = np.array(layer_alpha).transpose()
        self.save_alpha_data(self.save_alpha)
