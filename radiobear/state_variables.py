# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import print_function, absolute_import, division
import six


def init_state_variables(mode, **kwargs):
    state_vars = {'batch_mode': False,
                  'write_output_files': True,
                  'write_log_file': True,
                  'plot': None,  # Provided as a courtesy and for backward compatibility
                  'plot_atm': True,
                  'plot_bright': True,
                  'verbose': True,  # 0/None/False, 'normal'/True, 'loud'
                  'save_alpha': False,
                  'use_existing_alpha': False,
                  'scale_existing_alpha': False,
                  'normalize_weighting': True,
                  'output_type': 'frequency',  # or 'wavelength'
                  'log_directory': 'Logs',
                  'output_directory': 'Output',
                  'scratch_directory': 'Scratch'
                  }

    if mode == 'batch':
        state_vars['batch_mode'] = True
        state_vars['plot_atm'] = False
        state_vars['plot_bright'] = False
        state_vars['verbose'] = False
        state_vars['write_log_file'] = False
    elif mode == 'mcmc':
        state_vars['plot_atm'] = False
        state_vars['plot_bright'] = False
        state_vars['verbose'] = False
        state_vars['write_log_file'] = False
        state_vars['write_output_files'] = False
        state_vars['scale_existing_alpha'] = True
    elif mode == 'use_alpha':
        state_vars['plot_atm'] = False
        state_vars['plot_bright'] = False
        state_vars['verbose'] = True
        state_vars['use_existing_alpha'] = True
    elif mode == 'scale_alpha':
        state_vars['plot_atm'] = False
        state_vars['plot_bright'] = False
        state_vars['verbose'] = True
        state_vars['scale_existing_alpha'] = True

    if 'plot' in kwargs.keys():
        state_vars['plot_atm'] = kwargs['plot']
        state_vars['plot_bright'] = kwargs['plot']

    # update based on provided kwargs
    for k, v in six.iteritems(kwargs):
        if k in state_vars.keys():
            state_vars[k] = v
        else:
            print("'{}' keyword not found.".format(k))
            raise ValueError("Aborting since you probably wanted this keyword")

    # check various constraints
    only_one_allowed = ['batch_mode', 'save_alpha', 'use_existing_alpha', 'scale_existing_alpha']
    xxx = [state_vars[x] for x in only_one_allowed]
    if xxx.count(True) > 1:
        print(only_one_allowed)
        raise ValueError("Only one is allowed to be True")

    return state_vars


def set_state(state_class, set_mode='set', **kwargs):
    """
    set_mode:  'set' or 'init', if set, checks list
    """
    for k, v in six.iteritems(kwargs):
        if k in state_class.state_vars:
            setattr(state_class, k, v)
            if set_mode == 'set':
                print('Setting {} to {}'.format(k, v))
        else:
            if set_mode == 'set':
                print('state_var [{}] not found.'.format(k))
    if set_mode == 'init' and state_class.verbose == 'loud':
        show_state(state_class)


def show_state(state_class):
    print("{} state variables".format(str(state_class).strip('<').split()[0]))
    for k in state_class.state_vars:
        print('\t{}:  {}'.format(k, getattr(state_class, k)))
