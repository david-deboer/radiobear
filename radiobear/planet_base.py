#  This is the 'executive' class for planets
from __future__ import absolute_import, division, print_function
import numpy as np
import datetime
import os
import six
from argparse import Namespace
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


class PlanetBase:
    def __init__(self, name, mode='normal', config_file='config.par', **kwargs):
        """This is the base planet class to compute overall planetary emission.
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
        self.config_file = config_file

        print('Planetary modeling  (ver {})'.format(version.VERSION))
        if self.planet not in planet_list:
            print("{} not found.".format(self.planet))
            return

        # Set up state_variables
        self.mode = mode.lower()
        self.kwargs = state_variables.init_state_variables(self.mode, **kwargs)
        self.state_vars = self.kwargs.keys()
        state_variables.set_state(self, set_mode='init', **self.kwargs)

    def set_log(self):
        if self.write_log_file:
            runStart = datetime.datetime.now()
            logFile = '{}/{}_{}.log'.format(self.log_directory, self.planet, runStart.strftime("%Y%m%d_%H%M%S"))
            self.log = logging.LogIt(logFile)
            self.log.add(self.planet + ' start ' + str(runStart), self.verbose)
        else:
            self.log = None

    def set_data_return(self):
        self.data_return = data_handling.Data()
        self.data_return.set('log', self.log)

    def set_config(self):
        self.config_file = os.path.join(self.planet, self.config_file)
        if self.verbose:
            print('Reading config file:  ', self.config_file)
            print("\t'{}.config.display()' to see config parameters.".format(self.planet[0].lower()))
        self.config = config.planetConfig(self.planet, configFile=self.config_file, log=self.log)
        self.config.show()

    def state(self):
        state_variables.show_state(self)

    def init_atmos(self):
        #  ## Create atmosphere:  attributes are self.atmos.gas, self.atmos.cloud and self.atmos.layerProperty
        self.atmos = atmosphere.Atmosphere(self.planet, mode=self.mode, config=self.config, log=self.log, **self.kwargs)

    def init_alpha(self):
        #  ## Read in absorption modules:  to change absorption, edit files under /constituents'
        self.alpha = alpha.Alpha(mode=self.mode, config=self.config, log=self.log, **self.kwargs)
        self.alpha.reset_layers()

    def init_bright(self):
        #  ## Next compute radiometric properties - initialize bright and return data class
        self.bright = brightness.Brightness(mode=self.mode, log=self.log, **self.kwargs)

    def init_fIO(self):
        # ## Create fileIO class
        self.fIO = fileIO.FileIO(directory=self.output_directory)

    def set_bright_plots(self):
        if self.plot_bright:
            from radiobear.plotting import bright, data
            return bright.plots(self.bright), data.plots(self.data_return)
        else:
            return None

    def set_atm_plots(self):
        # ## Set plots
        if self.plot_atm:
            from radiobear.plotting import atm
            return atm.plots(self.atmos)
        else:
            return None

    def generate_freqs(self, freqs, freqUnit='GHz'):
        """
        Generates the frequencies to use.  Sets self.freqs and self.freqUnit.

        Parameters
        ----------
            freqs :
            freqUnit : str
                unit for freqs
        """
        if self.use_existing_alpha or self.scale_existing_alpha:
            freqs_read = np.load('{}/freqs.npy'.format(self.scratch_directory))
            freqs = [f for f in freqs_read]
            if self.verbose == 'loud':
                print("Setting frequencies to ", freqs)
            self.freqs = freqs
        else:
            self.freqs, self.freqUnit, reuse = self.set_freq(freqs, freqUnit)
            self.alpha.reset_layers()
        self.data_return.set('f', self.freqs)
        self.data_return.set('freqUnit', self.freqUnit)
        return reuse

    def generate_b(self, b=[0.0, 0.0], block=[1, 1]):
        """Runs the model to produce the brightness temperature, weighting functions etc etc
            freqs:  frequency request as set in set_freq.  If 'reuse' it won't recompute absorption/layer (allows many b)
            b:  "impact parameter" request as set in set_b
            block:  blocks to produce image (related to memory error...)"""
        #  ##Set b, etc
        self.b = self.set_b(b, block)
        self.block = block
        self.data_return.set('b', self.b)

    def set_image(self):
        block_postfix = '_'
        if self.data_type != 'image':
            return Namespace(true=False, block=block_postfix, imrow=[])

        if len(self.freqs) > 1:
            raise ValueError('Warning:  Image must be at only one frequency')
        if self.verbose == 'loud':
            print('imgSize = {} x {}'.format(self.imSize[0], self.imSize[1]))
        if abs(block[1]) > 1:
            block_postfix = '_{:02d}of{:02d}_'.format(self.block[0], abs(self.block[1]))
        return Namespace(true=True, block=block_postfix, imrow=[])

    def alpha_layers(self):
        self.alpha.get_layers(self.freqs, self.atmos)
        if self.save_alpha:
            self.alpha.write_alpha()

    def init_run(self):
        self.Tb = []
        self.rNorm = None
        self.tip = None
        self.rotate = None
        if self.verbose == 'loud':
            print('data_type = {}'.format(self.data_type))

    def get_bright(self, b, is_img):
        Tb = self.bright.single(self.freqs, self.atmos, b, self.alpha, self.config.orientation)
        if is_img.true:
            is_img.imrow.append(Tb[0])
            if not (i + 1) % self.imSize[0]:
                self.Tb.append(is_img.imrow)
                isimg.imrow = []
        else:
            self.Tb.append(Tb)
        if self.bright.travel is not None:
            if self.rNorm is None:
                self.rNorm = self.bright.travel.rNorm
            if self.tip is None:
                self.tip = self.bright.travel.tip
            if self.rotate is None:
                self.rotate = self.bright.travel.rotate
        return Tb

    def populate_data_return(self, runStart, runStop):
        self.data_return.set('start', runStart)
        self.data_return.set('stop', runStop)
        self.data_return.set('Tb', self.Tb)
        self.data_return.set('type', self.data_type)
        self.data_return.set('header', self.header)
        self.data_return.set('logfile', self.log.logfile)

    def set_header(self, missed_planet, run_start, run_stop):
        if missed_planet:
            self.header['orientation'] = '# orientation not set'
            self.header['aspect'] = '# aspect tip, rotate not set'
            self.header['rNorm'] = '# rNorm not set'
        else:
            self.header['orientation'] = '# orientation:   {}'.format(repr(self.config.orientation))
            self.header['aspect'] = '# aspect tip, rotate:  {:.4f}  {:.4f}'.format(utils.r2d(self.tip), utils.r2d(self.rotate))
            self.header['rNorm'] = '# rNorm: {}'.format(self.rNorm)
            if self.data_type == 'image':
                self.header['imgSize'] = '# imgSize: {}'.format(self.imSize)
                resolution = utils.r2asec(np.arctan(abs(self.b[1][0] - self.b[0][0]) * self.rNorm / self.config.distance))
                print('resolution = ', resolution)
                self.header['res'] = '# res:  {} arcsec'.format(resolution)
        self.header['data-type'] = '#* type:  {}'.format(self.data_type)
        self.header['gtype'] = '# gtype: {}'.format(self.config.gtype)
        self.header['radii'] = '# radii:  {:.1f}  {:.1f}  km'.format(self.config.Req, self.config.Rpol)
        self.header['distance'] = '# distance:  {} km'.format(self.config.distance)
        self.header['log-file:'] = '#* logfile: {}'.format(self.log.logfile)
        self.header['start'] = "#* start: {:%Y-%m-%d %H:%M:%S}".format(run_start)
        self.header['stop'] = "#* stop: {:%Y-%m-%d %H:%M:%S}".format(run_stop)

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
        reuse = False
        if len(self.freqs) == len(freqs):
            reuse = True
            for fslf, flcl in zip(sorted(self.freqs), sorted(freqs)):
                if (fslf - flcl) / fslf > 0.01:
                    reuse = False
                    break
        if reuse:
            freqs = self.freqs
            freqUnit = self.freqUnit
        else:
            if len(freqs) > 1:
                s = '{} in {} frequency steps ({} - {} {})'.format(self.planet, len(freqs), freqs[0], freqs[-1], utils.proc_unit(freqUnit))
            else:
                s = '{} at {} {}'.format(self.planet, freqs[0], utils.proc_unit(freqUnit))
            self.log.add(s, self.verbose)
            self.freqs = freqs
            self.freqUnit = utils.proc_unit(freqUnit)
        return freqs, self.freqUnit, reuse
