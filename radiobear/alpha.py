# -*- mode: python; coding: utf-8 -*-
"""Absorption calculations."""
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
import os
import sys
from argparse import Namespace
import numpy as np
from . import utils
from . import logging


class Alpha:
    """Read in absorption formalisms and compute layer absorption."""

    def __init__(self, idnum=0, config=None, log=None, load_formal=True, verbose=True, **kwargs):
        """
        Initialize class.

        Note that freq are all in GHz.
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
        self.constituentsAreAt = os.path.join(os.path.dirname(__file__), 'constituents')
        self.idnum = idnum
        self.reset_layers()
        self.alphafile = os.path.join(self.config.scratch_directory,
                                      'alpha{:04d}.npz'.format(self.idnum))
        if load_formal:
            self.setup_formalisms()
        self.saved_fields = ['ordered_constituents', 'alpha_data', 'freqs', 'P']
        self.memory = Namespace()

    def setup_formalisms(self):
        """Read in and setup formalism methods."""
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
        other_to_copy['co'] = ['coshape']
        self.other_dict = {}
        for absorber in self.ordered_constituents:
            self.other_dict[absorber] = {}
            if absorber in other_to_copy.keys():
                for oc in other_to_copy[absorber]:
                    try:
                        self.other_dict[absorber][oc] = getattr(self.config, oc)
                    except AttributeError:
                        print("The config for {} has no attribute {}".format(absorber, oc))

    def reset_layers(self):
        """Reset class parameters."""
        self.P = None
        self.freqs = None
        self.layers = None

    def save_alpha_data(self, save_type):
        """
        Write out the saved_fields.

        Parameters
        ----------
        save_type : str (other are ignored but don't error)
            'file' writes to self.alphafile
            'memory' writes to memory Namespace
        """
        self.tosave = np.array(self.tosave)
        if save_type == 'file':
            np.savez(self.alphafile,
                     ordered_constituents=self.ordered_constituents,
                     alpha_data=self.tosave,
                     freqs=self.freqs,
                     P=self.P)
        elif save_type == 'memory':
            self.memory.ordered_constituents = self.ordered_constituents
            self.memory.alpha_data = self.tosave
            self.memory.freqs = self.freqs
            self.memory.P = self.P

    def read_alpha_data(self, save_type):
        """
        Read the saved_fields into self.

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

    def total_layer_alpha(self, absorb, lscale):
        """
        Take the layer absorption profile and scale-sums it.

        Parameters
        ----------
        absorb : numpy array
                 layer absorptions
        lscale : number or dictionary
                 layer scale values
        """
        totalAbsorption = np.zeros_like(self.freqs)
        Nfreq = len(self.freqs)
        if self._save_alpha_memfil:
            save_absorb = np.zeros_like(absorb)
        if utils.isanynum(lscale):
            lscale = float(lscale)
            for i in range(Nfreq):
                totalAbsorption[i] = absorb[i].sum() * lscale
                if self._save_alpha_memfil:
                    for j in range(len(absorb[i])):
                        save_absorb[i, j] = absorb[i][j] * lscale
            if self._save_alpha_memfil:
                self.tosave.append(save_absorb)
                del save_absorb
            del absorb
            return totalAbsorption

        for i in range(Nfreq):
            for j in range(absorb.shape[1]):
                constituent = self.ordered_constituents[j]
                new_value = absorb[i, j]
                if constituent in lscale.keys():
                    new_value *= lscale[constituent]
                if self._save_alpha_memfil:
                    save_absorb[i][j] = new_value
                totalAbsorption[i] += new_value
        if self._save_alpha_memfil:
            self.tosave.append(save_absorb)
            del save_absorb
        del absorb
        return totalAbsorption

    def get_alpha_from_calc(self, freqs, T, P, gas, gas_dict, cloud, cloud_dict, units='invcm'):
        """
        Get the total absoprtion coefficient from gas.

        It assumes the correct frequency units, but maybe should correct that.
        """
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
        """
        Get the absorption coefficient for a layer.

        Either calculate from formalisms or read from saved data
        """
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
        """Get scale for a layer."""
        layer_scale = []
        if isinstance(scale, dict):
            for k, v in scale.items():
                if k not in self.ordered_constituents:
                    raise ValueError("{} not found as constituent for alpha".format(k))
                if len(v) != N:
                    raise ValueError("Incorrect scale for {}:  N {} vs {}".format(k, len(v), N))
            for lyr in range(N):
                layer_scale.append({})
                for k, v in scale.items():
                    layer_scale[lyr][k] = v[lyr]
        elif isinstance(scale, (list, np.ndarray)):
            if len(scale) != N:
                raise ValueError("Incorrect number of scale layers.")
            layer_scale = scale
        else:
            if utils.isanynum(scale):
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
        # Set up stuff
        self.reset_layers()
        self.alpha_data = []
        self.tosave = []
        self.freqs = freqs
        self.P = atm.gas[atm.config.C['P']]
        self._get_alpha_memfil = get_alpha in ['memory', 'file']
        self._save_alpha_memfil = save_alpha in ['memory', 'file']

        # Run stuff
        if self._get_alpha_memfil:
            self.read_alpha_data(get_alpha)
        numLayers = len(atm.gas[0])
        au = utils.alphaUnit
        self.log.add('{} layers'.format(numLayers), self.verbose)
        lscale = self.get_layer_scale(scale, numLayers)
        layer_alpha = []
        for layer in range(numLayers):
            layer_alpha.append(self.get_single_layer(freqs, layer, atm, lscale[layer], units=au))
        self.layers = np.array(layer_alpha).transpose()

        # Save/del stuff
        if self._save_alpha_memfil:
            self.save_alpha_data(save_alpha)
        del lscale, layer_alpha, self.alpha_data, self.tosave
