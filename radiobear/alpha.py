# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import os
import sys
import numpy as np
from . import utils
from . import config as pcfg
from . import state_variables
from . import logging
import six


def initialize_scalefile(fn, atm_pressure, constituents, values):
    with open(fn, 'w') as fp:
        fp.write("#p {}\n".format(' '.join(constituents)))
        for p in atm_pressure:
            fp.write("{} {}\n".format(p, ' '.join([str(x) for x in values])))


def write_scalefile(fn, columns, values):
    with open(fn, 'w') as fp:
        fp.write("#{}\n".format(' '.join(columns)))
        for i in range(len(values[columns[0]])):
            s = ''
            for col in columns:
                s += '{} '.format(values[col][i])
            fp.write(s.strip() + '\n')


def read_scalefile(fn):
    with open(fn, 'r') as fp:
        for line in fp:
            if line[0] == '#':  # Must the first line and be there!
                col = line.strip('#').strip().split()
                columns = [x.lower() for x in col]
                values = {}
                for col in columns:
                    values[col] = []
            else:
                data = [float(x) for x in line.split()]
                for col, d in zip(columns, data):
                    values[col].append(d)
    return columns, values


class Alpha:
    def __init__(self, mode='normal', config=None, log=None, **kwargs):
        """Reads in absorption formalisms
           Note that they are all in GHz"""

        kwargs = state_variables.init_state_variables(mode, **kwargs)
        self.state_vars = kwargs.keys()
        state_variables.set_state(self, set_mode='init', **kwargs)
        self.log = logging.setup(log)
        self.freqs = None
        self.constituentsAreAt = os.path.join(os.path.dirname(__file__), 'Constituents')
        self.absorb_layer_save_data = []

        # get config
        if config is None or isinstance(config, six.string_types):
            config = pcfg.planetConfig(self.planet, configFile=config, log=self.log)
        self.config = config

        self.alpha_data = None
        if self.use_existing_alpha or self.scale_existing_alpha:
            self.existing_alpha_setup()
        else:
            self.formalisms()
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
        self.layers = None
        self.absorb_layer_save_data = []

    def get_layers(self, freqs, atm):
        self.freqs = freqs
        numLayers = len(atm.gas[0])
        layerAlp = []
        self.log.add('{} layers'.format(numLayers), self.verbose)
        for layer in range(numLayers):
            layerAlp.append(self.getAlpha(freqs, layer, atm, units=utils.alphaUnit))
            if self.verbose == 'loud':
                print('\r\tAbsorption in layer {}   '.format(layer + 1), end='')
        self.layers = np.array(layerAlp).transpose()

    def write_alpha(self):
        constfile = os.path.join(self.scratch_directory, 'constituents')
        absorbfile = os.path.join(self.scratch_directory, 'absorb')
        freqfile = os.path.join(self.scratch_directory, 'freqs')
        np.savez(constfile, alpha_dict=self.config.constituent_alpha,
                 alpha_sort=self.ordered_constituents)
        self.absorb_layer_save_data = np.array(self.absorb_layer_save_data)
        np.save(absorbfile, self.absorb_layer_save_data)
        np.save(freqfile, self.freqs)

    def existing_alpha_setup(self):
        if self.alpha_data is None:
            self.alpha_data = np.load('{}/absorb.npy'.format(self.scratch_directory))
            condata = np.load('{}/constituents.npz'.format(self.scratch_directory))
            self.ordered_constituents = condata['alpha_sort']
        if self.scale_existing_alpha:
            self.scale_constituent_columns, self.scale_constituent_values =\
                                            read_scalefile(self.config.scale_file_name)

    def write_scale(self, fn):
        write_scalefile(fn, self.scale_constituent_columns, self.scale_constituent_values)

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

    def getAlpha(self, freqs, layer, atm, units='invcm'):
        """This is a wrapper to get the absorption coefficient, either from
           calculating from formalisms or reading from file"""
        if self.freqs is None and self.save_alpha:
            self.freqs = freqs
        if self.use_existing_alpha or self.scale_existing_alpha:
            if len(self.alpha_data) != len(atm.gas[0]):
                raise ValueError("Absorption and atmosphere don't agree")
        if self.use_existing_alpha:
            return self.get_alpha_from_file(freqs, layer, units)
        elif self.scale_existing_alpha:
            for c in self.scale_constituent_values:
                if len(self.scale_constituent_values[c]) != len(atm.gas[0]):
                    raise ValueError("Scaling and atmosphere don't agree")
            return self.scale_alpha_from_file(freqs, layer, units)
        else:
            P = atm.gas[atm.config.C['P']][layer]
            T = atm.gas[atm.config.C['T']][layer]
            gas = atm.gas[:, layer]
            cloud = atm.cloud[:, layer]
            return self.get_alpha_from_calc(freqs, T, P, gas, atm.config.C,
                                            cloud, atm.config.Cl, units)

    def get_alpha_from_file(self, freqs, layer, units='invcm'):
        totalAbsorption = np.zeros_like(freqs)
        for i in range(len(freqs)):
            totalAbsorption[i] = self.alpha_data[layer, i, -1]
        return totalAbsorption

    def scale_alpha_from_file(self, freqs, layer, units='invcm'):
        totalAbsorption = np.zeros_like(freqs)
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
            self.data_layer(absorb, totalAbsorption)
        del absorb
        return totalAbsorption

    def data_layer(self, absorb, totalAbsorption):
        this_layer = []
        for i in range(len(absorb)):
            for a2 in absorb[i]:
                this_layer.append(a2)
            this_layer.append(totalAbsorption[i])
        self.absorb_layer_save_data.append(this_layer)

    def state(self):
        state_variables.show_state(self)
