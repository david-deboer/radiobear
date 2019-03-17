import matplotlib.pyplot as plt
import numpy as np


def plot_raypath_stuff(r=None, b=None, ray=None):
    if ray is None:
        plt.figure('observer')
        plt.plot(r * b[0], r * b[1], '.', color='k')
    else:
        plt.figure('raypath-r')
        plt.plot(ray.r4ds, ray.ds)
        plt.figure('raypath-P')
        plt.semilogy(ray.ds, ray.P4ds)
        plt.axis(ymin=ray.P4ds[-1], ymax=ray.P4ds[0])


class shape_plots:
    def __init__(self, pltcls):
        self.pltcls = pltcls

    def plotShapes(self, planet, r, lat=90.0, delta_lng=0.0, gtypes=['ellipse', 'sphere', 'reference', 'gravity'], latstep='default'):
        colors = ['k', 'r', 'g', 'b', 'y', 'm', 'c']
        plt.figure('Shapes')
        for i, gtype in enumerate(gtypes):
            self.saveShape = [['oui']]
            rmag = self.pltcls.calcShape(planet, r, pclat=90.0, delta_lng=delta_lng, gtype=gtype, latstep=latstep)
            del self.pltcls.saveShape[0]  # remove nonsense first term
            self.pltcls.saveShape = np.array(self.pltcls.saveShape)
            if gtype != 'gravity':
                self.pltcls.saveShape = np.flipud(self.pltcls.saveShape)
            _y = []
            _z = []
            for v in self.pltcls.saveShape:
                _y.append(v[_Y])
                _z.append(v[_Z])
            plt.plot(_z, _y, color=colors[self.colorCounter % len(colors)], label=gtype)
            self.colorCounter += 1
            plt.axis('image')
            plt.legend()


class planet_plots:
    def __init__(self, pltcls):
        self.pltcls = pltcls

    def plot_profiles(self, b):
        plt.figure("Profile")
        Tbtr = np.transpose(self.pltcls.Tb)
        for j in range(len(self.pltcls.freqs)):
            frqs = ('%.2f %s' % (self.pltcls.freqs[j], self.pltcls.freqUnit))
            plt.plot(b, Tbtr[j], label=frqs)
        plt.legend()
        plt.xlabel('b')
        plt.ylabel('$T_B$ [K]')


class bright_plots:
    def __init__(self, pltcls):
        self.pltcls = pltcls

    def plot_W(self, freqs, integrated_W, normW4plot):
        # ####-----Weigthing functions
        plt.figure('INT_W')
        plt.plot(freqs, integrated_W)
        plt.title('Integrated weighting function')
        plt.xlabel('Frequency [GHz]')
        plt.figure('radtran')
        for i, f in enumerate(freqs):
            if normW4plot:
                wplot = self.pltcls.W[i] / np.max(self.pltcls.W[i])
            else:
                wplot = self.pltcls.W[i]
            if self.pltcls.output_type == 'frequency':
                label = (r'{:.1f} GHz').format(f)
            else:
                label = (r'{:.1f} cm').format(30.0 / f)
            plt.semilogy(wplot, self.pltcls.P, label=label, linewidth=3)
        plt.legend()
        plt.axis(ymin=100.0 * np.ceil(np.max(self.pltcls.P) / 100.0), ymax=1.0E-7 * np.ceil(np.min(self.pltcls.P) / 1E-7))
        plt.ylabel('P [bars]')

    def plot_Alpha(self, freqs):
        # ####-----Alpha
        plt.figure('alpha')
        for i, f in enumerate(freqs):
            if self.pltcls.output_type == 'frequency':
                label = (r'$\alpha$: {:.1f} GHz').format(f)
            else:
                label = (r'{:.1f} cm').format(30.0 / f)
            pl = list(self.pltcls.layerAlpha[i])
            del pl[0]
            plt.loglog(pl, self.pltcls.P, label=label)
        plt.legend()
        v = list(plt.axis())
        v[2] = 100.0 * np.ceil(np.max(self.pltcls.P) / 100.0)
        v[3] = 1.0E-7 * np.ceil(np.min(self.pltcls.P) / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')

    def plot_Tb(self, freqs):
        # ####-----Brightness temperature
        plt.figure('brightness')
        lt = '-'
        if (len(self.pltcls.Tb) == 1):
            lt = 'o'
        plt.plot(freqs, self.pltcls.Tb, lt)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Brightness temperature [K]')


class atm_plots:
    def __init__(self, pltcls):
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
