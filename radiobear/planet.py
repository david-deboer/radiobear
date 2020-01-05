# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import datetime
import os.path
from . import planet_base
from . import utils


class Planet(planet_base.PlanetBase):
    """
    This is the 'executive' function class to compute overall planetary emission.

    Parameters
    ----------
        name : str
            Planet name.  One of [Jupiter, Saturn, Uranus, Neptune]
        config_file : str
            Config file name.  If 'planet' sets to <name>/config.par
        run_atm : bool
            Flag to generate the atmosphere on setup
        load_formal : bool
            Flag
        verbose : str/bool
            Set verbosity
        setup : list
            List of modules to setup (excludes config, handled separately)
        kwargs
            E.g. 'plot_atm', etc (and other config parameters)
    """
    def __init__(self, name, config_file='config.par', run_atm=True,
                 load_formal=True, verbose=True,
                 setup=['log', 'data_return', 'atm', 'alpha', 'bright', 'fIO'],
                 **kwargs):
        self.load_formal = load_formal
        self.verbose = verbose
        super(Planet, self).__init__(name=name, config_file=config_file)

        # initialize and setup up modules/etc
        self.setup_config(**kwargs)
        for par in setup:
            getattr(self, 'setup_{}'.format(par))()

        self.log.add(self.planet, False)
        self.log.add(config_file, False)
        pars = self.config.show()
        self.log.add(pars, False)

        # run atmosphere
        if run_atm:
            self.atm_run(atm_run_type=self.config.atm_run_type)
        self.alpha_options = {'f': 'file', 'm': 'memory', 'n': 'none', 'c': 'none'}

    def run(self, freqs, b='disc', scale=False, get_alpha='calc',
            save_alpha='none', freqUnit='GHz', block=[1, 1], reuse_override='check'):
        """
        Runs the model to produce the brightness temperature, weighting functions etc etc

        Parameters
        ----------
        get_alpha : str
            String for what to do for alpha:  'file', 'memory', 'calc'
            If 'none', it calculates absorption
        save_alpha : str
            String for how to save alpha:  'file', 'memory', 'none'
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
        get_alpha = self.alpha_options[get_alpha[0].lower()]
        save_alpha = self.alpha_options[save_alpha[0].lower()]
        reuse = self.check_reuse(freqs, scale, get_alpha, reuse_override=reuse_override.lower())
        if not reuse:
            self.set_freqs(freqs=freqs, freqUnit=freqUnit)
            self.get_alpha = get_alpha
            self.scale = scale
        self.set_b(b=b, block=block)
        brtplt, datplt = self.set_bright_plots()
        is_img = self.set_image()
        C_timer = datetime.datetime.now()

        if not reuse:
            self.alpha_layers(freqs=self.freqs, atmos=self.atmos,
                              scale=scale, get_alpha=get_alpha, save_alpha=save_alpha)
            D_timer = datetime.datetime.now()
            print("Absoprtion calc took {:.1f} s".format(utils.timer(D_timer - C_timer)))

        #  Loop over b values
        self.init_run()
        runStart = datetime.datetime.now()
        self.log.add('Run start ' + str(runStart), False)
        for i, bv in enumerate(self.b):
            # Figure out which alpha to use for this b.  For now only one.
            j = self.map_b_to_atm(bv)
            is_img.i = i
            if self.verbose == 'loud':
                print('{} of {} (view {})  '.format(i + 1, len(self.b), bv), end='')
            self.bright_run(b=bv, freqs=self.freqs, atm=self.atmos[j], alpha=self.alpha[j],
                            is_img=is_img, brtplt=brtplt)
        runStop = datetime.datetime.now()
        missed_planet = self.rNorm is None
        self.set_header(missed_planet, runStart, runStop)
        self.log.add('Run stop ' + str(runStop), False)
        self.populate_data_return(runStart, runStop)
        print("RT calc took {:.1f} s".format(utils.timer(runStop - runStart)))

        #  ##Write output files
        if self.config.write_output_files:
            output_file = '{}_{}{}{}.dat'.format(self.planet, self.data_type, is_img.block,
                                                 runStart.strftime("%Y%m%d_%H%M%S"))
            output_file = os.path.join(self.config.output_directory, output_file)
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
