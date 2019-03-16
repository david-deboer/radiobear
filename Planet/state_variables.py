from __future__ import print_function, absolute_import, division


def init_state_variables(mode, **kwargs):
    state_vars = {'batch_mode': False,
                  'initialize': True,
                  'write_output_files': True,
                  'write_log_file': True,
                  'plot': True,
                  'verbose': True,  # 0/None/False, 'normal'/True, 'loud'
                  'generate_alpha': False,
                  'use_existing_alpha': False,
                  'scale_existing_alpha': False,
                  'output_type': 'frequency'  # or 'wavelength'
                  }

    if mode == 'batch':
        state_vars['batch_mode'] = True
        state_vars['self.plot'] = False
        state_vars['verbose'] = False
        state_vars['write_log_file'] = False
    elif mode == 'mcmc':
        state_vars['plot'] = False
        state_vars['verbose'] = False
        state_vars['write_log_file'] = False
        state_vars['write_output_files'] = False
        state_vars['scale_existing_alpha'] = True
    elif mode == 'use_alpha':
        state_vars['plot'] = False
        state_vars['verbose'] = True
        state_vars['use_existing_alpha'] = True
    elif mode == 'scale_alpha':
        state_vars['plot'] = False
        state_vars['verbose'] = True
        state_vars['scale_existing_alpha'] = True

    # update based on provided kwargs
    for k in kwargs:
        if k in state_vars.keys():
            state_vars[k] = kwargs[k]
        else:
            print("'{}' keyword not found.".format(k))
            raise ValueError("Aborting since you probably wanted this keyword")

    # check various constraints
    only_one_allowed = ['batch_mode', 'generate_alpha', 'use_existing_alpha', 'scale_existing_alpha']
    xxx = [state_vars[x] for x in only_one_allowed]
    if xxx.count(True) > 1:
        print(only_one_allowed)
        raise ValueError("Only one is allowed to be True")

    return state_vars
