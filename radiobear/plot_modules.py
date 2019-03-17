import matplotlib.pyplot as plt


class atm_plots:
    def _init__(self, pltcls):
        self.pltcls = pltcls

    def plot_diff(self):
        color_seq = ['b', 'k', 'r', 'm', 'c']
        clr = {}
        plt.figure('Scale difference')
        for i, gas in enumerate(col):
            if gas.lower() != 'p':
                clr[gas] = color_seq[i]
                present, g = self.pltcls.is_present(self.pltcls.gas[self.pltcls.config.C[gas.upper()]])
                plt.loglog(g, self.pltcls.gas[self.pltcls.config.C['P']], color=clr[gas], linestyle='--', label=gas + ' before')
        for i, gas in enumerate(col):
            if gas.lower() != 'p':
                present, g = self.pltcls.is_present(self.pltcls.gas[self.pltcls.config.C[gas.upper()]])
                plt.loglog(g, self.pltcls.gas[self.pltcls.config.C['P']], color=clr[gas], linestyle='-', label=gas + ' after')
        self.frame_plot('Fractional Abundance')

    def frame_plot(self, xlabel, show_legend=True):
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * np.ceil(self.pltcls.gas[self.pltcls.config.C['P']][-1] / 100.0)
        v[3] = 1.0E-7 * np.ceil(self.pltcls.gas[self.pltcls.config.C['P']][0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(xlabel)
        if show_legend:
            plt.legend()

    def plotTP(self, plot='auto'):
        """Plot the T-P profile"""
        if plot == 'auto':
            plt.figure(self.pltcls.planet + ': T-P')
        plt.title(self.pltcls.planet + ':  T-P profile')
        plt.loglog(self.pltcls.gas[self.pltcls.config.C['T']], self.pltcls.gas[self.pltcls.config.C['P']], label='T')
        self.frame_plot('T [K]', show_legend=False)

    def plotCloud(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the clouds"""
        if plot == 'auto':
            plt.figure(self.pltcls.planet + ': clouds')
        plt.title(self.pltcls.planet + ': clouds')
        for cloud in self.pltcls.config.Cl:
            present, cl = self.pltcls.is_present(self.pltcls.cloud[self.pltcls.config.Cl[cloud]])
            if cloud in dontPlot or not present:
                continue
            plt.loglog(cl, self.pltcls.cloud[self.pltcls.config.Cl['P']], label=cloud)
        self.frame_plot(r'Density [g/cm$^3$]')

    def plotGas(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
        """Plots the constituents"""
        if plot == 'auto':
            plt.figure(self.pltcls.planet + ': gas')
        plt.title(self.pltcls.planet + ': gas')
        for gas in self.pltcls.config.C:
            present, g = self.pltcls.is_present(self.pltcls.gas[self.pltcls.config.C[gas]])
            if gas in dontPlot or not present:
                continue
            plt.loglog(g, self.pltcls.gas[self.pltcls.config.C['P']], label=gas)
        self.frame_plot('Fractional Abundance')

    def plotProp(self, dontPlot=['Z', 'P', 'T'], plot='auto'):
        if plot == 'auto':
            plt.figure(self.pltcls.planet + ': Properties')
        plt.title(self.pltcls.planet + ': other')
        for other in self.pltcls.config.LP:
            present, g = self.pltcls.is_present(self.pltcls.layerProperty[self.pltcls.config.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g, self.pltcls.gas[self.pltcls.config.C['P']], label=other)
        self.frame_plot('Property Value')
