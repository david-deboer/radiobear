# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

from __future__ import print_function, absolute_import, division
import six


def init_state_variables(state_class, mode, **kwargs):
    state_dict = {'batch_mode': False,
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
        state_dict['batch_mode'] = True
        state_dict['plot_atm'] = False
        state_dict['plot_bright'] = False
        state_dict['verbose'] = False
        state_dict['write_log_file'] = False
    elif mode == 'mcmc':
        state_dict['plot_atm'] = False
        state_dict['plot_bright'] = False
        state_dict['verbose'] = False
        state_dict['write_log_file'] = False
        state_dict['write_output_files'] = False
        state_dict['scale_existing_alpha'] = True
    elif mode == 'use_alpha':
        state_dict['plot_atm'] = False
        state_dict['plot_bright'] = False
        state_dict['verbose'] = True
        state_dict['use_existing_alpha'] = True
    elif mode == 'scale_alpha':
        state_dict['plot_atm'] = False
        state_dict['plot_bright'] = False
        state_dict['verbose'] = True
        state_dict['scale_existing_alpha'] = True

    if 'plot' in kwargs.keys():
        state_dict['plot_atm'] = kwargs['plot']
        state_dict['plot_bright'] = kwargs['plot']

    # update based on provided kwargs
    for k, v in six.iteritems(kwargs):
        if k in state_dict.keys():
            state_dict[k] = v
        else:
            raise ValueError("Aborting since you probably wanted keyword {}".format(k))

    # check various constraints
    only_one_allowed = ['batch_mode', 'save_alpha', 'use_existing_alpha', 'scale_existing_alpha']
    xxx = [state_dict[x] for x in only_one_allowed]
    if xxx.count(True) > 1:
        print(only_one_allowed)
        raise ValueError("Only one is allowed to be True:  {}".format(','.join(only_one_allowed)))

    state_class.state_vars__status__ = 'init'
    state_class.state_vars = list(state_dict.keys())
    set_state(state_class, **state_dict)
    state_class.state_vars__status__ = 'set'


def set_state(state_class, **kwargs):
    """
    Make dictionary entries state_class attributes.
    """
    for k, v in six.iteritems(kwargs):
        if k in state_class.state_vars:
            setattr(state_class, k, v)
            if state_class.state_vars__status__ == 'set':
                print("Setting state variable {} to:  {}".format(k, v))
        else:
            print('state_var {} not found.'.format(k))

    if state_class.state_vars__status__ == 'init' and state_class.verbose == 'loud':
        show_state(state_class)


def show_state(state_class):
    print("{} state variables".format(str(state_class).strip('<').split()[0]))
    for k in state_class.state_vars:
        print('\t{}:  {}'.format(k, getattr(state_class, k)))
