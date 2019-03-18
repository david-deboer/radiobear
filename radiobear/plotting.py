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


def planet_profile(data):
    plt.figure('Profile')
    for i, f in enumerate(data.f):
        bvec = []
        Tvec = []
        for j, b in enumerate(data.b):
            bvec.append(np.sqrt(b[0] * b[0] + b[1] * b[1]))
            Tvec.append(data.Tb[j][i])
        plt.plot(bvec, Tvec, label=str(f))
    plt.legend()


def plot_raypath_stuff(b, ray, req=None, rpol=None):
    plt.figure('raypath-r')
    plt.plot(ray.r4ds, ray.ds)
    plt.xlabel('Radius [km]')
    plt.ylabel('Step [km]')
    plt.figure('raypath-P')
    plt.semilogy(ray.ds, ray.P4ds)
    plt.xlabel('Step [km]')
    plt.ylabel('Pressure [bar]')
    plt.axis(ymin=ray.P4ds[-1], ymax=ray.P4ds[0])
    plt.figure('Observer')
    if req is not None and rpol is not None:
        prod = req * rpol
        th = np.arange(0, (2.0 + 1.0 / 90.0) * np.pi, np.pi / 90.0)
        rrr = prod / np.sqrt(np.power(req * np.sin(th), 2.0) + np.power(rpol * np.cos(th), 2.0))
        x = rrr * np.cos(th)
        y = rrr * np.sin(th)
        plt.plot(x, y, 'k')
        plt.plot(req * b[0], req * b[1], 'ko')
    else:
        plt.plot(b[0], b[1], 'ko')
    plt.axis('image')


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


def plot_profiles(pltcls, b):
    plt.figure("Profile")
    Tbtr = np.transpose(pltcls.Tb)
    for j in range(len(pltcls.freqs)):
        frqs = ('%.2f %s' % (pltcls.freqs[j], pltcls.freqUnit))
        plt.plot(b, Tbtr[j], label=frqs)
    plt.legend()
    plt.xlabel('b')
    plt.ylabel('T_B [K]')


def frame_plot(P, xlabel, show_legend=True):
    v = list(plt.axis())
    if v[0] < 1E-10:
        v[0] = 1E-10
    v[2] = 100.0 * np.ceil(P[-1] / 100.0)
    v[3] = 1.0E-7 * np.ceil(P[0] / 1E-7)
    plt.axis(v)
    plt.ylabel('P [bars]')
    plt.xlabel(xlabel)
    if show_legend:
        plt.legend()


def plot_intW(freqs, int_W):
    # ####-----Weigthing functions
    plt.figure('INT_W')
    lt = '-'
    if len(freqs) == 1:
        lt = 'o'
    plt.plot(freqs, int_W, lt)
    plt.title('Integrated weighting function')
    plt.xlabel('Frequency [GHz]')


def plot_W(freqs, bright, normW4plot):
    plt.figure('radtran')
    for i, f in enumerate(freqs):
        norm = 1.0
        if normW4plot:
            norm = np.max(bright.W[i])
        wplot = bright.W[i] / norm
        if bright.output_type == 'frequency':
            label = ('{:.1f} GHz').format(f)
        else:
            label = ('{:.1f} cm').format(30.0 / f)
        plt.semilogy(wplot, bright.P, label=label, linewidth=3)
    frame_plot(bright.P, 'W', True)


def plot_Alpha(freqs, bright):
    # ####-----Alpha
    plt.figure('alpha')
    for i, f in enumerate(freqs):
        if bright.output_type == 'frequency':
            label = ('a: {:.1f} GHz').format(f)
        else:
            label = ('{:.1f} cm').format(30.0 / f)
        plt.loglog(bright.layerAlpha[i][1:], bright.P, label=label)
    frame_plot(bright.P, 'dB/km', True)


def plot_Tb(freqs, Tb):
    # ####-----Brightness temperature
    plt.figure('brightness')
    lt = '-'
    if (len(Tb) == 1):
        lt = 'o'
    plt.plot(freqs, Tb, lt)
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Brightness temperature [K]')


class atm_plots:
    def __init__(self, pltcls):
        self.pltcls = pltcls

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
