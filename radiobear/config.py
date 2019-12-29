# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import absolute_import, division, print_function
from . import utils
from . import logging
import json
import six
import os.path


def set_single_val(val, unit=None,
                   special={'true': True, 'false': False, 'none': None, 'null': None}):
    if isinstance(val, six.string_types):
        if val.lower() in special.keys():
            val = special[str(val).lower()]
        elif ',' in val:
            val = val.split(',')
        elif '.' in val or 'e' in val or 'E' in val:
            try:
                val = float(val)
            except ValueError:
                pass
        else:
            try:
                val = int(val)
            except ValueError:
                pass
    if isinstance(val, float):
        val = utils.convert_unit(val, unit)

    return val


default_config_file = os.path.join(os.path.dirname(__file__), 'default_config.json')


class planetConfig:
    def __init__(self, planet, configFile, log=None, printHelp=False):
        """reads in config file"""
        planet = planet.capitalize()
        self.planet = planet
        self.filename = configFile
        self.path = planet
        self.log = logging.setup(log)

        with open(default_config_file, 'r') as f:
            config_data = json.load(f)
        self.toks = config_data['toks']

        # Set defaults
        for tok in self.toks:
            val = self.toks[tok]['default'][self.planet]
            if isinstance(val, (six.string_types, float)):
                val = set_single_val(val, self.toks[tok]['unit'])
            setattr(self, self.toks[tok]['name'], val)
        self.setConfig(configFile)
        pars = self.show()
        self.log.add(planet, False)
        self.log.add(configFile, False)
        self.log.add(pars, False)

    def setConfig(self, configFile):
        """Reads in config files and updates after default set in __init__.
        These are all shown in display()"""
        if configFile is None:
            print('Config file not provided.  Using defaults.')
            return 0
        try:
            fp = open(configFile, 'r')
        except IOError:
            print(configFile, ' not found.  Using defaults.')
            return 0

        for line in fp:
            if line[0] in utils.commentChars or len(line) < 4:
                continue
            if '#' in line:
                line = line[:line.index('#')]
            data = line.split()
            tok = data[0].lower()
            del(data[0])
            if tok not in self.toks.keys():
                print('token {} not found'.format(tok))
                continue
            if isinstance(self.toks[tok]['default'][self.planet], (six.string_types, float, int)):
                unit = 'none'
                if len(data) == 2:
                    unit = data[1]
                if len(data):
                    val = set_single_val(data[0], unit)
                else:
                    print("{} didn't have an associated argument in config file.".format(tok))
                    print("Using default:  {}".format(self.toks[tok]['default'][self.planet]))
                    continue
            elif isinstance(self.toks[tok]['default'][self.planet], list):
                if len(data) == 1 and ',' in data[0]:
                    data = data[0].split(',')
                val = [set_single_val(x) for x in data]
            elif isinstance(self.toks[tok]['default'][self.planet], dict):
                val = {}
                for i, v in enumerate(data):
                    if ':' in v:
                        frml = v.split(':')
                        val[frml[0].lower().strip()] = set_single_val(frml[1].lower().strip())
                    else:
                        val[v.strip()] = i
            else:
                print("Incorrect type:  {} <{}>".
                      format(tok, type(self.toks[tok]['default'][self.planet])))
                continue
            setattr(self, self.toks[tok]['name'], val)
        fp.close()

        if 'DZ' not in self.C.keys():
            self.C['DZ'] = len(self.C.keys())
        if 'DZ' not in self.Cl.keys():
            self.Cl['DZ'] = len(self.Cl.keys())

        try:
            fp = open(self.zonal, 'r')
            self.vwlat = []
            self.vwdat = []
            for line in fp:
                data = line.split()
                self.vwlat.append(float(data[0]))
                self.vwdat.append(float(data[1]))
            fp.close()
        except IOError:
            self.vwlat = [0.0, 90.0]
            self.vwdat = [0.0, 0.0]

    def update_config(self, key, value, override=False):
        if not isinstance(key, list):
            key = [key]
        if not isinstance(value, list):
            value = [value]
        if len(key) != len(value):
            print("key/value pairs not matched.")
            return None
        for i, k in enumerate(key):
            if not override and k not in self.toks.keys():
                print("{} not in config keys and override is False")
                continue
            setattr(self, k, value[i])
            self.toks[k] = {}
            self.toks[k]['name'] = k
            self.toks[k]['unit'] = None
            self.toks[k]['help'] = None
            self.toks[k]['default'] = {'Jupiter': None, 'Saturn': None,
                                       'Uranus': None, 'Neptune': None}
        return (key, value)

    def show(self):
        """Returns string containing configuration"""
        s = 'Run parameters:\n'
        keys = list(self.toks.keys())
        keys.sort()
        for key in keys:
            s += '\t{:20s}:  {}\n'.format(key, str(getattr(self, self.toks[key]['name'])))
        return s
