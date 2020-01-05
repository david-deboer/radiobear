# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import absolute_import, division, print_function
from . import utils
import json
import os.path


def set_single_val(val, unit=None,
                   special={'true': True, 'false': False, 'none': None, 'null': None}):
    if isinstance(val, str):
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


default_config_files = [os.path.join(os.path.dirname(__file__), 'default_config.json'),
                        os.path.join(os.path.dirname(__file__), 'default_state.json')]


class planetConfig:
    """
    Initializes, reads in, updates and displays the configuration.
    """
    def __init__(self, planet, configFile):
        """reads in config file"""
        planet = planet.capitalize()
        self.planet = planet
        self.filename = configFile
        self.path = planet

        # Set defaults
        self.toks = {}
        for json_file in default_config_files:
            with open(json_file, 'r') as f:
                config_data = json.load(f)
            for tok, fval in config_data['toks'].items():
                self.toks[tok] = fval['name']
                pval = fval['default'][self.planet]
                if isinstance(pval, (str, float)):
                    pval = set_single_val(pval, fval['unit'])
                setattr(self, fval['name'], pval)
        self.setConfig(configFile)

    def __str__(self):
        return self.show()

    def setConfig(self, configFile):
        """Reads in config files and updates after defaults set in __init__.
        These are all shown in print(<instance>)"""
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
            try:
                preval = getattr(self, self.toks[tok])
            except ValueError:
                print('token {} ({}) not found'.format(tok, self.toks[tok]))
                continue
            if isinstance(preval, (str, float, int)):
                unit = 'none'
                if len(data) == 2:
                    unit = data[1]
                if len(data):
                    val = set_single_val(data[0], unit)
                else:
                    print("{} didn't have an associated argument in config file.".format(tok))
                    print("Retaining  {}".format(preval))
                    continue
            elif isinstance(preval, list):
                if len(data) == 1 and ',' in data[0]:
                    data = data[0].split(',')
                val = [set_single_val(x) for x in data]
            elif isinstance(preval, dict):
                val = {}
                for i, v in enumerate(data):
                    if ':' in v:
                        frml = v.split(':')
                        val[frml[0].lower().strip()] = set_single_val(frml[1].lower().strip())
                    else:
                        val[v.strip()] = i
            else:
                print("Incorrect type:  {} <{}>".
                      format(tok, type(preval)))
                continue
            setattr(self, self.toks[tok], val)
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

    def update_config(self, key=None, value=None, **kwargs):
        """
        Update the config file in a variety of ways.

        Parameters
        ----------
        key : any
            If dict: updates based on key/value
            If list: updates assuming value is same length list
            If str: updates based on key, value
        value : any
            See above for key
        override : bool
            Allow for new config entry to be made.
        kwargs : key = val
        """
        proc_key = []
        proc_val = []
        if key is None:
            pass
        elif isinstance(key, dict):
            for k, v in key.items():
                proc_key.append(k)
                proc_val.append(v)
        elif not isinstance(key, list):
            proc_key.append(key)
            proc_val.append(value)
        else:
            if len(key) != len(value):
                print("key/value pairs not matched.")
            else:
                for k, v in zip(key, value):
                    proc_key.append(k)
                    proc_val.append(v)
        for k, v in kwargs.items():
            proc_key.append(k)
            proc_val.append(v)

        for k, v in zip(proc_key, proc_val):
            if k in self.toks.keys():
                tok = self.toks[k]
            else:
                tok = k
            try:
                preval = getattr(self, tok)  # noqa
            except ValueError:
                print("{} not in config keys and override is False".format(tok))
                continue
            setattr(self, tok, set_single_val(v))

    def show(self):
        """Returns string containing configuration"""
        s = 'Run parameters:\n'
        keys = list(self.toks.keys())
        keys.sort()
        for key in keys:
            s += '\t{:20s}:  {}\n'.format(key, str(getattr(self, self.toks[key])))
        return s
