# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import datetime
import os.path
import six
from . import planet_base
from . import utils


class Planet(planet_base.PlanetBase):
    """
    This is the 'executive' function class to compute overall planetary emission.
    For both mode and kwargs look at state_variables.py

    Parameters
    ----------
        name : str
            Planet name.  One of [Jupiter, Saturn, Uranus, Neptune]
        mode : str
            Sets up for various special modes '[normal]/batch/mcmc/scale_alpha/use_alpha'
        config_file : str
            Config file name.  If 'planet' sets to <name>/config.par
        i_set : list
            List of str of modules to run on setup
        initialize : list
            List of str of modules to initialize on setup
        run_atmos : bool
            Flag to generate the atmosphere on setup
        kwargs
            'verbose' and 'plot_atm', etc (and other state_vars - see show_state())
    """
    def __init__(self, name, mode='normal', config_file='config.par',
                 set_log='log', set_config='config', set_data='data_return',
                 init_atm='atm', init_alpha='alpha', init_bright='bright', init_IO='fIO',
                 run_atmos=True, **kwargs):
        super(Planet, self).__init__(name=name, mode=mode, config_file=config_file, **kwargs)
        # set
        if isinstance(set_log, six.string_types) and set_log not in utils.negative:
            getattr(self, 'set_{}'.format(set_log))()
        if isinstance(set_config, six.string_types) and set_config not in utils.negative:
            getattr(self, 'set_{}'.format(set_config))()
        if isinstance(set_data, six.string_types) and set_data not in utils.negative:
            getattr(self, 'set_{}'.format(set_data))()
        # initialize
        if isinstance(init_atm, six.string_types) and init_atm not in utils.negative:
            getattr(self, 'init_{}'.format(init_atm))()
        if isinstance(init_alpha, six.string_types) and init_alpha not in utils.negative:
            getattr(self, 'init_{}'.format(init_alpha))()
        if isinstance(init_bright, six.string_types) and init_bright not in utils.negative:
            getattr(self, 'init_{}'.format(init_bright))()
        if isinstance(init_IO, six.string_types) and init_IO not in utils.negative:
            getattr(self, 'init_{}'.format(init_IO))()
        # run atmosphere
        if run_atmos:
            self.atm_run()

    def run(self, freqs, b=[0.0, 0.0], freqUnit='GHz', block=[1, 1]):
        """
        Runs the model to produce the brightness temperature, weighting functions etc etc

        Parameters
        ----------
        freqs : *
            frequency request as set in set_freq.
        b : *
            "impact parameter" request as set in set_b
        freqUnit : str
            unit of freqs
        block :  list
            blocks to produce image (related to memory error...)

        Returns
        -------
        data_return object
        """
        reuse = self.set_freqs(freqs=freqs, freqUnit=freqUnit)
        self.set_b(b=b, block=block)

        brtplt, datplt = self.set_bright_plots()
        is_img = self.set_image()
        C_timer = datetime.datetime.now()

        # For now just one profile, but can extend...
        if not reuse:
            self.alpha_layers()
        D_timer = datetime.datetime.now()
        print("Absoprtion calc took {:.1f} s".format(utils.timer(D_timer - C_timer)))

        #  Loop over b values
        self.init_run()
        runStart = datetime.datetime.now()
        self.log.add('Run start ' + str(runStart), False)
        for i, bv in enumerate(self.b):
            # Figure out which alpha to use for this b.  For now only one.
            if self.verbose == 'loud':
                print('{} of {} (view {})  '.format(i + 1, len(self.b), bv), end='')
            self.bright_run(b=bv, is_img=is_img, brtplt=brtplt, ibv=i)
        runStop = datetime.datetime.now()
        missed_planet = self.rNorm is None
        self.set_header(missed_planet, runStart, runStop)
        self.log.add('Run stop ' + str(runStop), False)
        self.populate_data_return(runStart, runStop)
        print("RT calc took {:.1f} s".format(utils.timer(runStop - runStart)))

        #  ##Write output files
        if self.write_output_files:
            output_file = '{}_{}{}{}.dat'.format(self.planet, self.data_type, is_img.block,
                                                 runStart.strftime("%Y%m%d_%H%M%S"))
            output_file = os.path.join(self.output_directory, output_file)
            if self.verbose == 'loud':
                print('\nWriting {} data to {}'.format(self.data_type, output_file))
            self.fIO.write(output_file, self.data_return)
        if brtplt is not None:
            brtplt.alpha()
            if self.data_type == 'spectrum' or self.data_type == 'profile' and len(freqs) > 1:
                datplt.Tb()
            if self.data_type == 'profile':
                datplt.profile()
            datplt.show()

        return self.data_return
