import matplotlib.pyplot as plt
import numpy as np


def plotTB(fn=None, xaxis='Frequency', xlog=False, justFreq=False, directory='Output', distance=4377233696.68):
    """plots brightness temperature against frequency and disc location:
           fn = filename to read (but then ignores directory) | '?', '.' or None | integer [None]
           xaxis = 'f[requency]' | 'w[avelength' ['freq']
           xlog = True | False [False]
           justFreq = True | False [False]
           directory = subdirectory for data (not used if filename given) ['Output']
           distance = distance for angular size plot in km [4377233696 km for Neptune]"""

    filename, Tb, f, wavel, b, xlabel, ylabels = read(fn=fn, directory=directory)
    title = filename.split('/')

    # Frequency plot
    plt.figure('TB')
    for i in range(len(b)):
        if xaxis[0].lower() == 'f':
            plotx = f
        else:
            plotx = wavel
            xlabel = 'Wavelength [cm]'
        if xlog:
            plt.semilogx(plotx, Tb[i], label=ylabels[i])
        else:
            plt.plot(plotx, Tb[i], label=ylabels[i])
    plt.xlabel(xlabel)
    plt.ylabel('Brightness Temperature [K]')
    plt.title(title[-1])

    if justFreq:
        return len(f)

    # # b plot
    plt.figure('b')
    for i in range(len(f)):
        plt.plot(b[:, 0], Tb[:, i], label=str(f[i]))
        plt.plot(b[:, 0], Tb[:, i], 'o')
    plt.legend()
    plt.title(title[-1])
    plt.xlabel('km')
    plt.ylabel('Brightness Temperature [K]')

    # b plot vs angle
    angle = []
    for r in b[:, 0]:
        angle.append((r / distance) * (180.0 / np.pi) * 3600.0)
    plt.figure('b_vs_angle')
    for i in range(len(f)):
        plt.plot(angle, Tb[:, i], label=str(f[i]))
        plt.plot(angle, Tb[:, i], 'o')
    plt.legend()
    plt.title(title[-1])
    plt.xlabel('arcsec')
    plt.ylabel('Brightness Temperature [K]')

    return len(b)


def plotObs(fn, cols=[0, 1, 2], color='b', marker='o', delimiter=None, comline='!'):
    try:
        fp = open(fn, 'r')
    except IOError:
        print(fn, ' not found')
        return 0
    data = []
    for line in fp:
        if comline in line[0:5]:
            continue
        dline = line.split(delimiter)
        if len(dline) < max(cols):
            continue
        drow = []
        for c in cols:
            drow.append(float(dline[c]))
        data.append(drow)
    data = np.array(data)
    plt.figure('ObsData')
    plt.semilogx(data[:, 0], data[:, 1], color=color, marker=marker)
    plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], color=color, marker=marker, ls='none')


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
        plt.ylabel('T_B [K]')


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
                label = ('{:.1f} GHz').format(f)
            else:
                label = ('{:.1f} cm').format(30.0 / f)
            plt.semilogy(wplot, self.pltcls.P, label=label, linewidth=3)
        plt.legend()
        plt.axis(ymin=100.0 * np.ceil(np.max(self.pltcls.P) / 100.0), ymax=1.0E-7 * np.ceil(np.min(self.pltcls.P) / 1E-7))
        plt.ylabel('P [bars]')

    def plot_Alpha(self, freqs):
        # ####-----Alpha
        plt.figure('alpha')
        for i, f in enumerate(freqs):
            if self.pltcls.output_type == 'frequency':
                label = ('a: {:.1f} GHz').format(f)
            else:
                label = ('{:.1f} cm').format(30.0 / f)
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
        self.frame_plot('Density [g/cm^3]')

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
