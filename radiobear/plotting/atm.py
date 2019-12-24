import matplotlib.pyplot as plt
import numpy as np


# ##############################################################################################################
#                                           ATMOSPHERE
# ##############################################################################################################
class plots:
    def __init__(self, atm):
        self.atm = atm

    def show(self):
        plt.show()

    def _frame_plot(self, xlabel, show_legend=True):
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * np.ceil(self.atm.gas[self.atm.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * np.ceil(self.atm.gas[self.atm.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(xlabel)
        if show_legend:
            plt.legend()

    def TP(self, plot='auto'):
        """Plot the T-P profile"""
        if plot == 'auto':
            plt.figure(self.atm.planet + ': T-P')
        plt.title(self.atm.planet + ':  T-P profile')
        plt.loglog(self.atm.gas[self.atm.config.C['T']], self.atm.gas[self.atm.config.C['P']],
                   label='T')
        self._frame_plot('T [K]', show_legend=False)

    def Cloud(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the clouds"""
        if plot == 'auto':
            plt.figure(self.atm.planet + ': clouds')
        plt.title(self.atm.planet + ': clouds')
        for cloud in self.atm.config.Cl:
            present, cl = self.atm.is_present(self.atm.cloud[self.atm.config.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl, self.atm.cloud[self.atm.config.Cl['P']], label=cloud)
        self._frame_plot('Density [g/cm^3]')

    def Gas(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the constituents"""
        if plot == 'auto':
            plt.figure(self.atm.planet + ': gas')
        plt.title(self.atm.planet + ': gas')
        for gas in self.atm.config.C:
            present, g = self.atm.is_present(self.atm.gas[self.atm.config.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g, self.atm.gas[self.atm.config.C['P']], label=gas)
        self._frame_plot('Fractional Abundance')

    def Properties(self, dontPlot=['Z', 'P', 'T'], plot='auto'):
        if plot == 'auto':
            plt.figure(self.atm.planet + ': Properties')
        plt.title(self.atm.planet + ': other')
        for other in self.atm.config.LP:
            present, g = self.atm.is_present(self.atm.property[self.atm.config.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g, self.atm.gas[self.atm.config.C['P']], label=other)
        self._frame_plot('Property Value')
