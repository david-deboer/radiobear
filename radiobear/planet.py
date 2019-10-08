#  This is the 'executive' class for planets
from __future__ import absolute_import, division, print_function
import numpy as np
import datetime
import os
import six
from . import atmosphere
from . import config
from . import alpha
from . import brightness
from . import data_handling
from . import utils
from . import fileIO
from . import state_variables
from . import logging
from . import version
from . import planet_base


class Planet(planet_base.PlanetBase):
    def __init__(self, name, mode='normal', config_file='config.par', **kwargs):
        """This is the 'executive' function class to compute overall planetary emission.
           For both mode and kwargs look at state_variables.py
           Inputs:
                name:  'Jupiter', 'Saturn', 'Uranus', 'Neptune'
                config_file:  config file name.  If 'planet' sets to <name>/config.par
                mode:  sets up for various special modes '[normal]/batch/mcmc/scale_alpha/use_alpha'
                kwargs: 'verbose' and 'plot_amt', etc (and other state_vars - see show_state())"""
        super(Planet, self).__init__(name=name, mode=mode, config_file=config_file, **kwargs)
        self.set_logfile()
        self.set_configfile()
        self.set_data_return()
        if self.initialize:
            self.init_atm()
            self.init_alpha()
            self.init_brightness()
            self.init_IO()

    def run(self, freqs='reuse', b=[0.0, 0.0], freqUnit='GHz', block=[1, 1]):
        """Runs the model to produce the brightness temperature, weighting functions etc etc
            freqs:  frequency request as set in set_freq.  If 'reuse' it won't recompute absorption/layer (allows many b)
            b:  "impact parameter" request as set in set_b
            freqUnit:  unit that freqs is in
            block:  blocks to produce image (related to memory error...)"""

        self.atmos.run()
        atmplt = self.atm_plots()
        if atmplot is not None:
            atmp


        self.generate_freqs(freqs=freqs, freqUnit=freqUnit)
        self.generate_b(b=b, block=block)

        brtplt = self.bright_plots()

        if self.data_type == 'image' and len(self.freqs) > 1:
            raise ValueError('Warning:  Image must be at only one frequency')
        if self.verbose == 'loud':
            print('data_type = {}'.format(self.data_type))

        # ##Initialize other stuff
        self.Tb = []
        btmp = ''
        self.rNorm = None
        self.tip = None
        self.rotate = None
        if self.data_type == 'image':  # We now treat it as an image at one frequency
            if self.verbose == 'loud':
                print('imgSize = {} x {}'.format(self.imSize[0], self.imSize[1]))
            imtmp = []
            if abs(block[1]) > 1:
                btmp = '_{:02d}of{:02d}'.format(block[0], abs(block[1]))

        #  ##Start b loop
        runStart = datetime.datetime.now()
        self.log.add('Run start ' + str(runStart), self.verbose)
        self.data_return.set('start', runStart)
        for i, bv in enumerate(b):
            if self.verbose == 'loud':
                print('{} of {} (view {})  '.format(i + 1, len(b), bv), end='')
            Tbt = self.bright.single(freqs, self.atmos, bv, self.alpha, self.config.orientation)
            if self.bright.travel is not None:
                if self.rNorm is None:
                    self.rNorm = self.bright.travel.rNorm
                if self.tip is None:
                    self.tip = self.bright.travel.tip
                if self.rotate is None:
                    self.rotate = self.bright.travel.rotate
            if self.data_type == 'image':
                imtmp.append(Tbt[0])
                if not (i + 1) % self.imSize[0]:
                    self.Tb.append(imtmp)
                    imtmp = []
            else:
                self.Tb.append(Tbt)
            if self.plot_bright:
                brtplt.raypath()
                brtplt.observer(b=bv, req=self.config.Req, rpol=self.config.Rpol)
                brtplt.intW()
                brtplt.W(self.normalize_weighting)
        self.data_return.set('Tb', self.Tb)
        missed_planet = self.rNorm is None

        if self.generate_alpha:
            self.alpha.complete_generate_alpha()

        runStop = datetime.datetime.now()
        self.log.add('Run stop ' + str(runStop), self.verbose)
        self.data_return.set('stop', runStop)
        self.data_return.set('type', self.data_type)
        self.set_header(missed_planet)
        self.data_return.set('header', self.header)
        self.data_return.set('logfile', self.log.logfile)

        #  ##Write output files
        if self.write_output_files:
            output_file = '{}/{}_{}{}_{}.dat'.format(self.output_directory, self.planet, self.data_type, btmp, runStart.strftime("%Y%m%d_%H%M%S"))
            if self.verbose == 'loud':
                print('\nWriting {} data to {}'.format(self.data_type, output_file))
            self.fIO.write(output_file, self.data_return)
        if self.plot_bright:
            brtplt.alpha()
            datplt = data.plots(self.data_return)
            if self.data_type == 'spectrum' or self.data_type == 'profile' and len(freqs) > 1:
                datplt.Tb()
            if self.data_type == 'profile':
                datplt.profile()
            datplt.show()

        return self.data_return
