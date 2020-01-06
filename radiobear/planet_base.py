# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
import numpy as np
import datetime
import os
import sys
from argparse import Namespace
from . import utils
from . import set_utils
from . import version


class PlanetBase:
    """
    This is the base planet class to compute overall planetary emission.

    Parameters
    ----------
        name : str
            One of [Jupiter, Saturn, Uranus, Neptune]
        config_file : str
            Config file name.  Uses <planet>/<config_file>
    """
    planet_list = ['Jupiter', 'Saturn', 'Neptune', 'Uranus']

    def __init__(self, name, config_file='config.par'):
        self.planet = name.capitalize()
        self.config_file = config_file

        self.header = {}
        self.freqs = []
        self.freqUnit = None
        self.b = None
        self.data_type = None
        self.imSize = None
        self.bmap_loaded = False

        self.version = version.VERSION
        self.version_notes = version.version_notes[version.VERSION]
        print('Planetary modeling  (ver {})'.format(self.version))
        if self.planet not in self.planet_list:
            print("{} not found.".format(self.planet))
            return

    def setup_config(self, **kwargs):
        """
        Instantiates and reads the config file
        """
        from . import config as pcfg
        self.config_file = os.path.join(self.planet, self.config_file)
        if self.verbose:
            print('Reading config file:  ', self.config_file)
            print("\t'print({}.config)' to see config parameters."
                  .format(self.planet[0].lower()))
        self.config = pcfg.planetConfig(self.planet, configFile=self.config_file)
        self.config.update_config(**kwargs)
        self.config.show()
        sys.path.insert(0, self.config.path)

    def setup_log(self):
        """
        Sets up self.log if self.write_log_file is True
        """
        from . import logging
        if self.config.write_log_file:
            runStart = datetime.datetime.now()
            logFile = '{}/{}_{}.log'.format(self.config.log_directory, self.planet,
                                            runStart.strftime("%Y%m%d_%H%M%S"))
            self.log = logging.LogIt(logFile)
            self.log.add(self.planet + ' start ' + str(runStart), self.verbose)
        else:
            self.log = None

    def setup_data_return(self):
        """
        Instantiates and logs the data_return class
        """
        from . import data_handling
        self.data_return = data_handling.Data()
        self.data_return.set('log', self.log)

    def setup_atm(self):
        """
        Instantiates atmosphere.  Attributes are:
            self.atmos[].gas, self.atmos[].cloud and self.atmos[].property
        """
        from . import atmosphere
        if not isinstance(self.config.gasFile, list):
            self.config.gasFile = [self.config.gasFile]
        if not isinstance(self.config.cloudFile, list):
            self.config.cloudFile = [self.config.cloudFile]
        N = len(self.config.gasFile)
        self.atmos = []
        for i in range(N):
            self.atmos.append(atmosphere.Atmosphere(self.planet, idnum=i,
                              config=self.config, log=self.log, verbose=self.verbose))

    def setup_alpha(self):
        """
        Instantiates absorption modules.  To change absorption, edit files under /constituents'
        """
        from . import alpha
        N = len(self.atmos)
        self.alpha = []
        for i in range(N):
            self.alpha.append(alpha.Alpha(idnum=i, config=self.config, log=self.log,
                                          load_formal=self.load_formal, verbose=self.verbose))

    def setup_bright(self):
        """
        Instantiates brightness module.
        """
        from . import brightness
        self.bright = brightness.Brightness(config=self.config, log=self.log, verbose=self.verbose)

    def setup_fIO(self):
        """
        Instantiates fileIO class
        """
        from . import fileIO
        self.fIO = fileIO.FileIO(directory=self.config.output_directory)

    def set_bright_plots(self):
        """
        If config.plot_bright is True, reads in brightness plotting modules.
        """
        if self.config.plot_bright:
            from radiobear.plotting import bright, data
            return bright.plots(self.bright), data.plots(self.data_return)
        else:
            return None, None

    def set_atm_plots(self, atmos):
        """
        If config.plot_atm is True, reads in atmosphere plotting modules.
        """
        if self.config.plot_atm:
            from radiobear.plotting import atm
            return atm.plots(atmos)
        else:
            return None

    def set_freqs(self, freqs, freqUnit='GHz'):
        """
        Sets the frequencies to use.  Sets self.freqs and self.freqUnit.

        Parameters
        ----------
            freqs : (see set_utils.set_freq)
            freqUnit : str
                unit for freqs
        """
        self.header['freqs'] = '# freqs request: {} {}'.format(str(freqs), freqUnit)
        freqs, freqUnit = set_utils.set_freq(freqs, freqUnit)
        self.data_return.set('f', freqs)
        self.data_return.set('freqUnit', freqUnit)
        return freqs, freqUnit

    def set_b(self, b=[0.0, 0.0], block=[1, 1]):
        """
        Sets the "impact parameter".  Sets self.b, self.block, self.data_type, self.imSize

        Parameter
        ---------
            b : (see set_utils.set_b)
                "impact parameter" request as set in set_b
            block : list
                blocks to produce image (related to memory error...)
        """
        self.header['b'] = '# b request:  {}  {}'.format(str(b), str(block))
        rv = set_utils.set_b(b, block, Rpol=self.config.Rpol, Req=self.config.Req)
        self.b = rv.b
        self.block = rv.block
        self.data_type = rv.data_type
        self.imSize = rv.imSize
        self.data_return.set('b', self.b)

    def map_b_to_atm(self, b=None):
        """
        Given the b index and value, returns the appropriate atmosphere index
        """
        if self.config.bmapmodule is None or b is None:
            return 0
        if not self.bmap_loaded:
            __import__(self.config.bmapmodule)
            bmapModule = sys.modules[self.config.bmapmodule]
            self.bmap_loaded = True
        return bmapModule.bmap(b=b)

    def set_image(self):
        """
        Returns the Namespace pertaining to data_type == image
        """
        block_postfix = '_'
        if self.data_type != 'image':
            return Namespace(true=False, block=block_postfix, imrow=[])

        if len(self.freqs) > 1:
            raise ValueError('Warning:  Image must be at only one frequency')
        if self.verbose == 'loud':
            print('imgSize = {} x {}'.format(self.imSize[0], self.imSize[1]))
        if abs(self.block[1]) > 1:
            block_postfix = '_{:02d}of{:02d}_'.format(self.block[0], abs(self.block[1]))
        return Namespace(true=True, i=None, block=block_postfix, imrow=[])

    def alpha_layers(self, freqs, atmos, scale=False, get_alpha='calc', save_alpha='none'):
        """
        Computes the layer absorption for all atmospheres.  If save_alpha is set,
        it will write the profiles to file or memory.
        """
        for i, atm in enumerate(atmos):
            if scale:
                if isinstance(scale, bool) and scale:
                    print("read self.scale_by and act")
                elif isinstance(scale, str):
                    print("read in file for scale")
                elif isinstance(scale, dict) or isinstance(scale, float):
                    pass
                else:
                    scale = False
            else:
                scale = False
            self.alpha[i].reset_layers()
            self.alpha[i].get_layers(freqs=freqs,
                                     atm=atm,
                                     scale=scale,
                                     get_alpha=get_alpha,
                                     save_alpha=save_alpha)

    def init_run(self):
        """
        Initializes class variables for computing brightness temperature
        """
        self.Tb = []
        self.rNorm = None
        self.tip = None
        self.rotate = None
        if self.verbose == 'loud':
            print('data_type = {}'.format(self.data_type))

    def atm_run(self, atm_run_type='std'):
        for atm in self.atmos:
            getattr(atm, atm_run_type)()
            atmplt = self.set_atm_plots(atmos=atm)
            if atmplt is not None:
                atmplt.TP()
                atmplt.Gas()
                atmplt.Cloud()
                atmplt.Properties()
                atmplt.show()

    def bright_run(self, b, freqs, atm, alpha, is_img, brtplt):
        """
        Computes the brightness temperature for that "b" and updates self.Tb.

        Parameters
        ----------
        b : list [float<1.0, float<1.0]
            Current "impact parameter"
        freq : list
            List of frequencies
        atm : Atmosphere class
            Atmosphere to use
        alpha : Alpha class
            Absorption to use
        is_img : Namespace
            Contains parameters if image.
        brtplt : class
            Content for plotting

        Returns
        -------
        float
            Brightness temperature at that b
        """
        Tb = self.bright.single(b, freqs, atm, alpha, self.config.orientation)
        if is_img.true:
            is_img.imrow.append(Tb[0])
            if not (is_img.i + 1) % self.imSize[0]:
                self.Tb.append(is_img.imrow)
                is_img.imrow = []
        else:
            self.Tb.append(Tb)
        if self.bright.travel is not None:
            if self.rNorm is None:
                self.rNorm = self.bright.travel.rNorm
            if self.tip is None:
                self.tip = self.bright.travel.tip
            if self.rotate is None:
                self.rotate = self.bright.travel.rotate
        if brtplt is not None:
            brtplt.raypath()
            brtplt.observer(b=b, req=self.config.Req, rpol=self.config.Rpol)
            brtplt.intW()
            brtplt.W(self.config.normalize_weighting)

        return Tb

    def populate_data_return(self, run_start, run_stop):
        """
        Populates the data_return instance with the run data.

        Parameters
        ----------
        run_start : datetime
            Time run started
        run_stop : datetime
            Time run ended
        """
        self.data_return.set('start', run_start)
        self.data_return.set('stop', run_stop)
        self.data_return.set('Tb', self.Tb)
        self.data_return.set('type', self.data_type)
        self.data_return.set('header', self.header)
        self.data_return.set('logfile', self.log.logfile)

    def check_reuse(self, freqs, scale, get_alpha, save_alpha, reuse_override='check'):
        if reuse_override == 'true':
            return True
        if reuse_override == 'false':
            return False
        # Check get_alpha
        if get_alpha != self.get_alpha:
            return False
        if save_alpha != self.save_alpha:
            return False
        # Check freqs
        if len(freqs) != len(self.freqs):
            return False
        for fslf, flcl in zip(sorted(self.freqs), sorted(freqs)):
            if (fslf - flcl) / fslf > 0.01:
                return False
        # Check scale
        if not isinstance(scale, type(self.scale)):
            return False
        if type(scale) in (float, int):
            if scale != self.scale:
                return False
        if isinstance(scale, (list, np.ndarray)):
            if len(scale) != len(self.scale):
                return False
            for a, b in zip(scale, self.scale):
                if (a - b) / a > 0.001:
                    return False
        if isinstance(scale, dict):
            alist = sorted(list(scale.keys()))
            blist = sorted(list(self.scale.keys()))
            if len(alist) != len(blist):
                return False
            for a, b in zip(alist, blist):
                if a != b:
                    return False
            for k in alist:
                if len(scale[k]) != len(self.scale[k]):
                    return False
                for a, b in zip(scale[k], self.scale[k]):
                    if (a - b) / a > 0.001:
                        return False
        return True

    def set_header(self, missed_planet, run_start, run_stop):
        """
        Populates the header dictionary.

        Parameters
        ----------
        missed_planet : bool
            Flag whether ray missed the planet or not
        run_start : datetime
            Time run started
        run_stop : datetime
            Time run ended
        """
        if missed_planet:
            self.header['orientation'] = '# orientation not set'
            self.header['aspect'] = '# aspect tip, rotate not set'
            self.header['rNorm'] = '# rNorm not set'
        else:
            self.header['orientation'] = '# orientation:   {}'.format(repr(self.config.orientation))
            self.header['aspect'] = '# aspect tip, rotate:  {:.4f}  {:.4f}'.format(
                                                                            utils.r2d(self.tip),
                                                                            utils.r2d(self.rotate))
            self.header['rNorm'] = '# rNorm: {}'.format(self.rNorm)
            if self.data_type == 'image':
                self.header['imgSize'] = '# imgSize: {}'.format(self.imSize)
                resolution = utils.r2asec(np.arctan(abs(self.b[1][0] - self.b[0][0]) *
                                          self.rNorm / self.config.distance))
                print('resolution = ', resolution)
                self.header['res'] = '# res:  {} arcsec'.format(resolution)
        self.header['data-type'] = '#* type:  {}'.format(self.data_type)
        self.header['gtype'] = '# gtype: {}'.format(self.config.gtype)
        self.header['radii'] = '# radii:  {:.1f}  {:.1f}  km'.format(self.config.Req,
                                                                     self.config.Rpol)
        self.header['distance'] = '# distance:  {} km'.format(self.config.distance)
        self.header['log-file:'] = '#* logfile: {}'.format(self.log.logfile)
        self.header['start'] = "#* start: {:%Y-%m-%d %H:%M:%S}".format(run_start)
        self.header['stop'] = "#* stop: {:%Y-%m-%d %H:%M:%S}".format(run_stop)
