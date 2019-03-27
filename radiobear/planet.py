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
        self.data_type = None
        self.imSize = None

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
            logFile = '{}/{}_{}.log'.format(self.log_directory, self.planet, runStart.strftime("%Y%m%d_%H%M%S"))
            self.log = logging.LogIt(logFile)
            self.log.add(self.planet + ' start ' + str(runStart), self.verbose)
        else:
            self.log = None
        self.data_return = data_handling.Data()
        self.data_return.set('log', self.log)

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

        # ## Create fileIO class
        self.fIO = fileIO.FileIO(directory=self.output_directory)

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
        if self.data_type == 'image' and len(freqs) > 1:
            print('Warning:  Image must be at only one frequency')
            print('Using {} {}'.format(freqs[0], freqUnit))
            self.freqs = list(freqs[0])
            freqs = self.freqs
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

    def set_header(self, missed_planet):
        if missed_planet:
            self.header['orientation'] = '# orientation not set'
            self.header['aspect'] = '# aspect tip, rotate not set'
            self.header['rNorm'] = '# rNorm not set'
        else:
            self.header['orientation'] = '# orientation:   {}'.format(repr(self.config.orientation))
            self.header['aspect'] = '# aspect tip, rotate:  {:.4f}  {:.4f}'.format(utils.r2d(self.tip), utils.r2d(self.rotate))
            self.header['rNorm'] = '# rNorm: {}'.format(self.rNorm)
            self.header['data_type'] = '#* type:  {}'.format(self.data_type)
            if self.data_type == 'image':
                self.header['imgSize'] = '# imgSize: {}'.format(self.imSize)
                resolution = utils.r2asec(np.arctan(abs(self.b[1][0] - self.b[0][0]) * self.rNorm / self.config.distance))
                print('resolution = ', resolution)
                self.header['res'] = '# res:  {} arcsec'.format(resolution)
        self.header['gtype'] = '# gtype: {}'.format(self.config.gtype)
        self.header['radii'] = '# radii:  {:.1f}  {:.1f}  km'.format(self.config.Req, self.config.Rpol)
        self.header['distance'] = '# distance:  {} km'.format(self.config.distance)
        self.header['log-file:'] = '#* logfile: {}'.format(self.log.logfile)
        self.header['start'] = "#* start: {:%Y-%m-%d %H:%M:%S}".format(self.data_return.start)
        self.header['stop'] = "#* stop: {:%Y-%m-%d %H:%M:%S}".format(self.data_return.stop)

    def set_b(self, b, block):
        """Process b request.
           Sets data_type ('image', 'spectrum', 'profile'), and imSize

           Returns list of b-pairs

           Parameters:
           ------------
           b:   list of  pairs -> returns that list
                one pair -> returns that as a list
                float -> returns a full image grid
                'disc' (or 'disk') -> returns [[0.0, 0.0]]
                'stamp:bres:xmin,xmax,ymin,ymax' -> returns grid of postage stamp
                'start:stop:step[<angle] -> string defining line
                'n1,n2,n3[<angle]' -> csv list of b magnitudes
           block:  image block as pair, e.g. [4, 10] is "block 4 of 10"
           """
        self.header['b'] = '# b request:  {}  {}'.format(str(b), str(block))
        # Deal with strings
        if isinstance(b, six.string_types):
            b = b.lower()
            if b.startswith('dis'):
                self.data_type = 'spectrum'
                return [b]
            if b.startswith('stamp'):
                bres = float(b.split(':')[1])
                bext = [float(x) for x in b.split(':')[2].split(',')]
                b = []
                for x in np.arange(bext[0], bext[1] + bres / 2.0, bres):
                    for y in np.arange(bext[2], bext[3] + bres / 2.0, bres):
                        b.append([y, x])
                xbr = len(np.arange(bext[2], bext[3] + bres / 2.0, bres))
                self.imSize = [xbr, len(b) / xbr]
                self.data_type = 'image'
            else:
                b = b.split('<')
                angle_b = 0.0 if len(b) == 1 else utils.d2r(float(b[1]))
                if ',' in b[0]:
                    mag_b = [float(x) for x in b[0].split(',')]
                elif ':' in b[0]:
                    mag_b = [float(x) for x in b[0].split(':')]
                    mag_b = np.arange(mag_b[0], mag_b[1] + mag_b[2] / 2.0, mag_b[2])
                ab = self.config.Rpol / self.config.Req
                rab = ab / np.sqrt(np.power(np.sin(angle_b), 2.0) + np.power(ab * np.cos(angle_b), 2.0))
                b = []
                for v in mag_b:
                    if v < 0.99 * rab:
                        b.append([v * np.cos(angle_b), v * np.sin(angle_b)])
                self.data_type = 'profile'
            return b
        if isinstance(b, float):  # this generates a grid at that spacing and blocking
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
            self.data_type = 'image'
            return b
        shape_b = np.shape(b)
        if len(shape_b) == 1:
            self.data_type = 'spectrum'
            return [b]
        else:
            self.data_type = 'spectrum' if shape_b[0] < 5 else 'profile'
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
        self.header['freqs'] = '# freqs request: {} {}'.format(str(freqs), freqUnit)
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
