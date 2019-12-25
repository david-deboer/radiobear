# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
import numpy as np
import scipy.special as ss
import os.path
from . import utils
from . import raypath as ray
from . import state_variables
from . import logging


class Brightness():

    def __init__(self, mode='normal', log=None, **kwargs):
        """This calculates the brightness temperature of the planets.
           It must be used with atmosphere and alpha"""
        state_variables.init_state_variables(self, mode, **kwargs)
        self.log = logging.setup(log)

    def state(self):
        state_variables.show_state(self)

    def single(self, b, freqs, atm, alpha, orientation=None, taulimit=20.0):
        """This computes the brightness temperature along one ray path"""

        print("BRIGHT28:  ",np.shape(alpha.layers), alpha.layers)
        disc_average = utils.b_type(b).startswith('dis')
        if disc_average:
            b = [0.0, 0.0]
        self.alpha = alpha
        self.freqs = freqs
        self.b = b

        if self.alpha.layers is None:
            self.alpha.get_layers(freqs, atm)
        # get path lengths (ds_layer) vs layer number (num_layer)
        #     currently frequency independent refractivity
        print_meta = (self.verbose == 'loud')
        travel = ray.compute_ds(atm, b, orientation, gtype=None, verbose=print_meta)
        self.travel = travel
        if travel.ds is None:
            print('Off planet')
            self.Tb = []
            for j in range(len(freqs)):
                self.Tb.append(utils.T_cmb)
            return self.Tb

        # set and initialize arrays
        integrated_W = [0.0 for f in freqs]
        self.tau = [[0.0 for f in freqs]]
        self.Tb_lyr = [[0.0 for f in freqs]]
        self.W = [[0.0 for f in freqs]]

        P_layers = atm.gas[atm.config.C['P']]
        T_layers = atm.gas[atm.config.C['T']]
        z_layers = atm.gas[atm.config.C['Z']]
        self.P = [P_layers[travel.layer4ds[0]]]
        self.z = [z_layers[travel.layer4ds[0]]]

        for i in range(len(travel.ds) - 1):
            ds = travel.ds[i] * utils.Units[utils.atmLayerUnit] / utils.Units['cm']
            taus = []
            Ws = []
            Tbs = []
            ii = travel.layer4ds[i]
            ii1 = travel.layer4ds[i + 1]
            T1 = T_layers[ii1]
            T0 = T_layers[ii]
            self.P.append((P_layers[ii] + P_layers[ii1]) / 2.0)
            self.z.append((z_layers[ii] + z_layers[ii1]) / 2.0)

            if self.alpha.layers is None:
                print("is None at ", i)
            for j, f in enumerate(freqs):
                if not alpha.config.Doppler:
                    a1 = self.alpha.layers[j][ii1]
                    a0 = self.alpha.layers[j][ii]
                else:
                    print("\n\nDoppler currently broken since the getAlpha call is different.")
                    fshifted = [[f / travel.doppler[i]], [f / travel.doppler[i + 1]]]
                    print('\rdoppler corrected frequency at layer', i, end='')
                    a1 = alpha.getAlpha(fshifted[0], T_layers[ii1], P_layers[ii1], atm.gas[:, ii1],
                                        atm.config.C, atm.cloud[:, ii1],
                                        atm.config.Cl, units=utils.alphaUnit)
                    a0 = alpha.getAlpha(fshifted[1], T_layers[ii], P_layers[ii], atm.gas[:, ii],
                                        atm.config.C, atm.cloud[:, ii],
                                        atm.config.Cl, units=utils.alphaUnit)
                dtau = (a0 + a1) * ds / 2.0
                taus.append(self.tau[i][j] + dtau)         # this is tau_(i+1)
                if disc_average:
                    Ws.append(2.0 * a1 * ss.expn(2, taus[j]))  # this is W_(i+1) for disc average
                else:
                    Ws.append(a1 * np.exp(-taus[j]))  # this is W_(i+1) for non disc average
                integrated_W[j] += (Ws[j] + self.W[i][j]) * ds / 2.0
                dTb = (T1 * Ws[j] + T0 * self.W[i][j]) * ds / 2.0
                Tbs.append(self.Tb_lyr[i][j] + dTb)
            self.tau.append(taus)
            self.W.append(Ws)
            self.Tb_lyr.append(Tbs)

        # final spectrum
        self.Tb = []
        for j in range(len(freqs)):
            top_Tb_lyr = self.Tb_lyr[-1][j]
            if top_Tb_lyr < utils.T_cmb:
                top_Tb_lyr = utils.T_cmb
            else:
                top_Tb_lyr /= integrated_W[j]  # Normalize by integrated weights (makes assumptions)
                if integrated_W[j] < 0.96 and self.verbose:
                    print("Weight correction at {:.2f} is {:.4f} (showing below 0.96)"
                          .format(freqs[j], integrated_W[j]))
            self.Tb.append(top_Tb_lyr)
        self.tau = np.array(self.tau).transpose()
        self.W = np.array(self.W).transpose()
        self.Tb_lyr = np.array(self.Tb_lyr).transpose()
        self.P = np.array(self.P)
        self.z = np.array(self.z)
        self.integrated_W = np.array(integrated_W)

        del taus, Tbs, Ws, travel
        return self.Tb

    def savertm(self, tag=None, path=None):
        if tag is None:
            filename = None
        else:
            filename = 'alpha_' + tag + '.out'
        self.saveAlpha(filename, self.output_directory)
        if tag is None:
            filename = None
        else:
            filename = 'wgt_' + tag + '.out'
        self.saveWeight(filename, path)
        if tag is None:
            filename = None
        else:
            filename = 'tau_' + tag + '.out'
        self.saveTau(filename, path)
        if tag is None:
            filename = None
        else:
            filename = 'tblayer_' + tag + '.out'
        self.saveTblayer(filename, path)

    def saveit(self):
        for i, f in enumerate(self.freqs):
            filename = 'pawtt_{:.3f}.out'.format(f)
            fp = open(filename, 'w')
            print("{}:  Pressure, alpha, weight, tau, Tb".format(filename))
            for j in range(len(self.P)):
                s = '{}\t'.format(repr(self.P[j]))
                s += '{}\t'.format(repr(self.alpha.layers[i][j]))
                s += '{}\t'.format(repr(self.W[i][j]))
                s += '{}\t'.format(repr(self.tau[i][j]))
                s += '{}\n'.format(repr(self.Tb_lyr[i][j]))
                fp.write(s)
            fp.close()

    def saveAlpha(self, filename=None, path='.'):
        if filename is None:
            filename = 'alpha.out'
        filename = os.path.join(path, filename)
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += '{:.2f}\t'.format(f)
        s += 'GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = ('{}\t{:.2f}\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += '{}\t'.format(repr(self.alpha.layers[i][j]))
            s += '\n'
            fp.write(s)
        s = ('{} ({} x {})').format(filename, i + 1, j + 1)

    def saveWeight(self, norm=False, filename=None, path='.'):
        if filename is None:
            filename = 'wgt.out'
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += ('{:.2f}\t').format(f)
        s = s.strip() + 'GHz\n'
        fp.write(s)
        scale = []
        for i in range(len(self.freqs)):
            if norm:
                scale.append(np.max(self.W[i]))
            else:
                scale.append(1.0)
        for j in range(len(self.P)):
            s = ('{}\t{:.2f}\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += ('{}\t').format(repr(self.W[i][j] / scale[i]))
            s = s.strip() + '\n'
            fp.write(s)
        s = ('{} ({} x {})').format(filename, i + 1, j + 1)
        return s

    def saveTau(self, filename=None, path='.'):
        if filename is None:
            filename = 'tau.out'
        os.path.join(path, filename)
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += '{:.2f}\t'.format(f)
        s += 'GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = ('{}\t{:.2f}\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += ('{}\t').format(repr(self.tau[i][j]))
            s += '\n'
            fp.write(s)
        s = ('{} ({} x {})').format(filename, i + 1, j + 1)
        return s

    def saveTblayer(self, filename=None, path='.'):
        if filename is None:
            filename = 'tblayer.out'
        os.path.join(path, filename)
        fp = open(filename, 'w')
        s = '#P  \tz  \t'
        for f in self.freqs:
            s += ('{:.2f}\t').format(f)
        s += 'GHz\n'
        fp.write(s)
        for j in range(len(self.P)):
            s = ('{}\t{:.2f}\t').format(repr(self.P[j]), self.z[j])
            for i in range(len(self.freqs)):
                s += ('{}\t').format(repr(self.Tb_lyr[i][j]))
            s += '\n'
            fp.write(s)
        s = ('{} ({} x {})').format(filename, i + 1, j + 1)
        return s
