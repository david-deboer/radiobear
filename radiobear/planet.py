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


class Planet:
    def __init__(self, name, mode='normal', config_file='config.par', **kwargs):
        """This is the 'executive function class to compute overall planetary emission.
           For both mode and kwargs look at state_variables.py
           Inputs:
                name:  'Jupiter', 'Saturn', 'Uranus', 'Neptune'
                config_file:  config file name.  If 'planet' sets to <name>/config.par
                mode:  sets up for various special modes '[normal]/batch/mcmc/scale_alpha/use_alpha'
                kwargs: 'verbose' and 'plot_amt', etc (and other state_vars - see show_state())"""

        planet_list = ['Jupiter', 'Saturn', 'Neptune', 'Uranus']
        self.planet = name.capitalize()
        self.header = {}
        self.freqs = None
        self.freqUnit = None
        self.b = None

        print('Planetary modeling  (ver {})'.format(version.VERSION))
        if self.planet not in planet_list:
            print("{} not found.".format(self.planet))
            return

        # Set up state_variables
        kwargs = state_variables.init_state_variables(mode.lower(), **kwargs)
        self.state_vars = kwargs.keys()
        state_variables.set_state(self, set_mode='init', **kwargs)
        self.mode = mode
        self.kwargs = kwargs

        #  ##Set up log file
        if self.write_log_file:
            runStart = datetime.datetime.now()
            logFile = '{}/{}_{}.log'.format(self.log_directory, self.planet, runStart.strftime("%Y%m%d_%H%M"))
            self.log = logging.LogIt(logFile)
            self.log.add(self.planet + ' start ' + str(runStart), self.verbose)
        else:
            self.log = None

        #  ## Get config
        config_file = os.path.join(self.planet, config_file)
        if self.verbose:
            print('Reading config file:  ', config_file)
            print("\t'{}.config.display()' to see config parameters.".format(name[0].lower()))
        self.config = config.planetConfig(self.planet, configFile=config_file, log=self.log)
        self.config.show()

        if self.initialize:
            self.initialize_run()

    def state(self):
        state_variables.show_state(self)

    def initialize_run(self):
        #  ## Create atmosphere:  attributes are self.atmos.gas, self.atmos.cloud and self.atmos.layerProperty
        self.atmos = atmosphere.Atmosphere(self.planet, mode=self.mode, config=self.config, log=self.log, **self.kwargs)
        self.atmos.run()

        #  ## Read in absorption modules:  to change absorption, edit files under /constituents'
        self.alpha = alpha.Alpha(mode=self.mode, config=self.config, log=self.log, **self.kwargs)

        #  ## Next compute radiometric properties - initialize bright and return data class
        self.bright = brightness.Brightness(mode=self.mode, log=self.log, **self.kwargs)
        self.data_return = data_handling.DataReturn()

        # ## Create fileIO class
        self.fIO = fileIO.FileIO(self.output_type, self.output_directory)

        # ## Set plots
        if self.plot_atm:
            from radiobear.plotting import atm
            atmplt = atm.plots(self.atmos)
            atmplt.TP()
            atmplt.Gas()
            atmplt.Cloud()
            atmplt.Properties()
            atmplt.show()

    def run(self, freqs='reuse', b=[0.0, 0.0], freqUnit='GHz', block=[1, 1]):
        """Runs the model to produce the brightness temperature, weighting functions etc etc
            freqs:  frequency request as set in set_freq.  If 'reuse' it won't recompute absorption/layer (allows many b)
            b:  "impact parameter" request as set in set_b
            freqUnit:  unit that freqs is in
            block:  blocks to produce image (related to memory error...)"""

        #  ##Set freqs
        if self.use_existing_alpha or self.scale_existing_alpha:
            freqs_read = np.load('{}/freqs.npy'.format(self.scratch_directory))
            freqs = [f for f in freqs_read]
            if self.verbose == 'loud':
                print("Setting frequencies to ", freqs)
        reuse = False
        if isinstance(freqs, six.string_types) and freqs.lower() == 'reuse':
            if self.freqs is None:
                raise ValueError('Must set frequencies.')
            reuse = True
            freqs = self.freqs
            freqUnit = self.freqUnit
        else:
            freqs, freqUnit = self.set_freq(freqs, freqUnit)
            self.bright.resetLayers()
        self.data_return.set('f', freqs)
        self.data_return.set('freqUnit', freqUnit)
        if self.plot_bright:
            from radiobear.plotting import bright, data
            brtplt = bright.plots(self.bright)

        #  ##Set b, etc
        b = self.set_b(b, block)
        self.data_return.set('b', b)
        if self.outType == 'image' and len(freqs) > 1:
            print('Warning:  Image must be at only one frequency')
            print('Using {} {}'.format(freqs[0], freqUnit))
            self.freqs = list(freqs[0])
            freqs = self.freqs
        if self.verbose == 'loud':
            print('outType = {}'.format(self.outType))

        # ##Initialize other stuff
        self.Tb = []
        btmp = ''
        self.rNorm = None
        self.tip = None
        self.rotate = None
        if self.outType == 'image':  # We now treat it as an image at one frequency
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
                print('{} of {} (view [{:.4f}, {:.4f}])  '.format(i + 1, len(b), bv[0], bv[1]), end='')
            Tbt = self.bright.single(freqs, self.atmos, bv, self.alpha, self.config.orientation, discAverage=(self.bType == 'disc'))
            if self.bright.travel is not None:
                if self.rNorm is None:
                    self.rNorm = self.bright.travel.rNorm
                if self.tip is None:
                    self.tip = self.bright.travel.tip
                if self.rotate is None:
                    self.rotate = self.bright.travel.rotate
            if self.outType == 'image':
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
        self.data_return.header = self.header
        missed_planet = self.rNorm is None

        if self.generate_alpha:
            self.alpha.complete_generate_alpha()

        #  ##Write output files
        if self.write_output_files:
            outputFile = '{}/{}_{}{}_{}.dat'.format(self.output_directory, self.planet, self.outType, btmp, runStart.strftime("%Y%m%d_%H%M"))
            if self.verbose == 'loud':
                print('\nWriting {} data to {}'.format(self.outType, datFile))
            self.set_header(missed_planet)
            self.fIO.write(outputFile, self.outType, freqs, freqUnit, b, self.Tb, self.header)
        if self.plot_bright:
            brtplt.alpha()
            datplt = data.plots(self.data_return)
            if self.outType == 'spectrum' or self.outType == 'profile' and len(freqs) > 1:
                datplt.Tb()
            if self.outType == 'profile':
                datplt.profile()
            datplt.show()

        runStop = datetime.datetime.now()
        self.log.add('Run stop ' + str(runStop), self.verbose)
        self.data_return.set('stop', runStop)
        return self.data_return

    def set_header(self, missed_planet):
        if missed_planet:
            self.header['res'] = '# res not set\n'
            self.header['orientation'] = '# orientation not set\n'
            self.header['aspect'] = '# aspect tip, rotate not set\n'
            self.header['rNorm'] = '# rNorm not set\n'
        else:
            self.header['orientation'] = '# orientation:   {}\n'.format(repr(self.config.orientation))
            self.header['aspect'] = '# aspect tip, rotate:  {:.4f}  {:.4f}\n'.format(utils.r2d(self.tip), utils.r2d(self.rotate))
            self.header['rNorm'] = '# rNorm: {}\n'.format(self.rNorm)
            if self.bType:
                self.header['bType'] = '# bType:  {}\n'.format(self.bType)
            if self.outType:
                self.header['outType'] = '# outType:  {}\n'.format(self.outType)
                if self.outType == 'image':
                    self.header['imgSize'] = '# imgSize: {}\n'.format(self.imSize)
                    resolution = utils.r2asec(np.arctan(abs(self.b[1][0] - self.b[0][0]) * self.rNorm / self.config.distance))
                    print('resolution = ', resolution)
                    self.header['res'] = '# res:  {} arcsec\n'.format(resolution)
        self.header['gtype'] = '# gtype: {}\n'.format(self.config.gtype)
        self.header['radii'] = '# radii:  {:.1f}  {:.1f}  km\n'.format(self.config.Req, self.config.Rpol)
        self.header['distance'] = '# distance:  {} km\n'.format(self.config.distance)

    def set_b(self, b, block):
        """Process b request.
           bType (see below), outType ('image', 'spectrum', 'profile'), and imSize get set as attributes
           b is list of ordered pairs on return

           Parameters:
           ------------
           b:   list of ordered pairs -> returns that list (bType='points')
                float -> returns a full grid (bType='image')
                'disc' (or 'disk') -> returns [[0.0, 0.0]] (bType='disc')
                'stamp' -> returns grid of postage stamp (bType='stamp')
                string defining line: csv list of values or range as <start>:<stop>:<step>
                                      optional ',angle=DEG' [defaults to 0.0]
                                      (bType='line')
           block:  image block as pair, e.g. [4, 10] is "block 4 of 10"
           """
        self.header['b'] = '# b request:  {}  {}\n'.format(str(b), str(block))
        self.imSize = None
        self.outType = 'spectrum'
        # Do some pre-processing to handle line vs single point and ndarrays
        if len(np.shape(b)) == 1 and len(b) > 2:
            b = ','.join([str(x) for x in b])
        if isinstance(b, six.string_types) and len(b.split(',')) == 2:
            b = b.split(',')

        # Handle lists/arrays
        if len(np.shape(b)) > 0:
            if len(np.shape(b)) == 1 and len(b) == 2:
                b = [b]
            if len(np.shape(b)) == 2:  # This is just a list of b pairs, so return.
                self.bType = 'points'
                return b
            else:
                raise ValueError("Invalid b request.")

        if isinstance(b, float):  # this generates a grid at that spacing and blocking
            self.bType = 'image'
            self.outType = 'image'
            grid = -1.0 * np.flipud(np.arange(b, 1.5 + b, b))
            grid = np.concatenate((grid, np.arange(0.0, 1.5 + b, b)))
            # get blocks
            bsplit = len(grid) / abs(block[1])
            lastRow = block[0] / abs(block[1])
            if abs(block[1]) == 1:
                lastRow = 0
            b = []
            for i in range(int(bsplit + lastRow)):
                ii = i + int((block[0] - 1) * bsplit)
                vrow = grid[ii]
                for vcol in grid:
                    b.append([vcol, vrow])
            self.imSize = [len(grid), len(b) / len(grid)]
            return b

        if not isinstance(b, six.string_types):
            raise ValueError("Invalid b request.")
        self.bType = b.lower()

        if self.bType == 'disc' or self.bType == 'disk':
            self.bType = 'disc'
            return [[0.0, 0.0]]

        if self.bType == 'stamp':
            self.outType = 'image'
            print('Setting to postage stamp.  Need more information')
            bres = float(raw_input('...Input postage stamp resolution in b-units:  '))
            bx = [float(x) for x in raw_input('...Input bx_min, bx_max:  ').split(',')]
            by = [float(x) for x in raw_input('...Input by_min, by_max:  ').split(',')]
            b = []
            for x in np.arange(bx[0], bx[1] + bres / 2.0, bres):
                for y in np.arange(by[0], by[1] + bres / 2.0, bres):
                    b.append([y, x])
            xbr = len(np.arange(by[0], by[1] + bres / 2.0, bres))
            self.imSize = [xbr, len(b) / xbr]
            return b

        # It is a line request
        line = {'mag_b': [], 'angle_b': 0.0, 'range': ':' in self.bType}
        cmd = self.bType.split(',')
        self.bType = 'line'
        self.outType = 'profile'
        for v in cmd:
            if '<' in v:
                line['angle_b'] = utils.d2r(float(v.split('<')[1]))
                v = v.split('<')[0]
            if ':' in v:
                start, stop, step = [float(x) for x in v.split(':')]
                line['mag_b'] = np.arange(start, stop + step / 2.0, step)
            if not line['range'] and len(v):
                line['mag_b'].append(float(v))
        b = []
        for v in line['mag_b']:
            ang = line['angle_b']
            ab = self.config.Rpol / self.config.Req
            rrr = ab / np.sqrt(np.power(np.sin(ang), 2.0) + np.power(ab * np.cos(ang), 2.0))
            if abs(v) < rrr:
                b.append([v * np.cos(ang), v * np.sin(ang)])
        return b

    def set_freq(self, freqs, freqUnit):
        """ Process frequency request.
            Return a list converted from freqUnit to processingFreqUnit and reassigns freqUnit procUnit.

            Parameters:
            ------------
            freqs:  list -> returns that list
                    csv list of values -> converts to list
                    float/int -> returns that one value as a list
                    '<start>:<stop>:<step>' -> returns arange of that
                    '<start>;<stop>;<nstep>' -> returns logspace of that
                    '<filename>' -> returns loadtxt of that
            freqUnit:  frequency unit of supplied freqs
        """
        self.header['freqs'] = '# freqs request: {} {}\n'.format(str(freqs), freqUnit)
        # ## Process frequency range "request"
        if isinstance(freqs, list):
            pass
        elif isinstance(freqs, np.ndarray):
            freqs = list(freqs)
        elif isinstance(freqs, six.string_types) and ',' in freqs:
            freqs = [float(x) for x in freqs.split(',')]
        elif isinstance(freqs, (float, int)):
            freqs = [float(freqs)]
        elif isinstance(freqs, six.string_types) and ':' in freqs:
            fstart, fstop, fstep = [float(x) for x in freqs.split(':')]
            freqs = list(np.arange(fstart, fstop + fstep / 2.0, fstep))
        elif isinstance(freqs, six.string_types) and ';' in freqs:
            fstart, fstop, nstep = [float(x) for x in freqs.split(';')]
            freqs = list(np.logspace(np.log10(fstart), np.log10(fstop), nstep))
        elif isinstance(freqs, six.string_types):
            freqs = list(np.loadtxt(freqs))
        else:
            raise ValueError('Invalid format for frequency request')

        for i in range(len(freqs)):
            freqs[i] = utils.convert_unit(freqs[i], freqUnit)
        if len(freqs) > 1:
            s = '{} in {} frequency steps ({} - {} {})'.format(self.planet, len(freqs), freqs[0], freqs[-1], utils.proc_unit(freqUnit))
        else:
            s = '{} at {} {}'.format(self.planet, freqs[0], utils.proc_unit(freqUnit))
        self.log.add(s, self.verbose)
        self.freqs = freqs
        self.freqUnit = utils.proc_unit(freqUnit)
        return freqs, self.freqUnit
