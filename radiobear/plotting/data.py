import matplotlib.pyplot as plt
import numpy as np


# ##############################################################################################################
#                                          GENERAL FILE PLOTTING
# ##############################################################################################################
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


# ##############################################################################################################
#                                              BRIGHTNESS
# ##############################################################################################################
class bright_plots:
    def __init__(self, bright):
        self.bright = bright

    def plot_Alpha(self):
        # ####-----Alpha
        plt.figure('alpha')
        for i, f in enumerate(self.bright.freqs):
            if self.bright.output_type == 'frequency':
                label = ('a: {:.1f} GHz').format(f)
            else:
                label = ('{:.1f} cm').format(30.0 / f)
            plt.loglog(self.bright.layerAlpha[i][1:], self.bright.P, label=label)
        self._frame_plot('dB/km', True)

    def plot_intW(self):
        # ####-----Weigthing functions
        plt.figure('INT_W')
        lt = '-'
        if len(self.bright.freqs) == 1:
            lt = 'o'
        plt.plot(self.bright.freqs, self.bright.integrated_W, lt)
        plt.title('Integrated weighting function')
        plt.xlabel('Frequency [GHz]')

    def plot_W(self, normW4plot):
        plt.figure('radtran')
        for i, f in enumerate(self.bright.freqs):
            norm = 1.0
            if normW4plot:
                norm = np.max(self.bright.W[i])
            wplot = self.bright.W[i] / norm
            if self.bright.output_type == 'frequency':
                label = ('{:.1f} GHz').format(f)
            else:
                label = ('{:.1f} cm').format(30.0 / f)
            plt.semilogy(wplot, self.bright.P, label=label, linewidth=3)
        self._frame_plot('W', True)

    def plot_raypath_stuff(self, b, req=None, rpol=None):
        ray = self.bright.travel
        plt.figure('raypath-r')
        if ray.ds is not None:
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

    def _frame_plot(self, xlabel, show_legend=True):
        v = list(plt.axis())
        if v[0] < 1E-10:
            v[0] = 1E-10
        v[2] = 100.0 * np.ceil(self.bright.P[-1] / 100.0)
        v[3] = 1.0E-7 * np.ceil(self.bright.P[0] / 1E-7)
        plt.axis(v)
        plt.ylabel('P [bars]')
        plt.xlabel(xlabel)
        if show_legend:
            plt.legend()


class data_plots:
    def __init__(self, data):
        self.data = data

    def planet_profile(self):
        plt.figure('Profile')
        b = self.data.b.transpose()
        bvec = np.sqrt(b[0] * b[0] + b[1] * b[1])
        Tvec = np.transpose(self.data.Tb)
        for i, f in enumerate(self.data.f):
            plt.plot(bvec, Tvec[i], label="{:.2f} {}".format(f, self.data.freqUnit))
        plt.legend()
        plt.xlabel('b')
        plt.ylabel('T_B [K]')

    def plot_Tb(self):
        # ####-----Brightness temperature
        plt.figure('brightness')
        for j, b in enumerate(self.data.b):
            Tplt = self.data.Tb[j]
            lt = '-'
            if (len(Tplt) == 1):
                lt = 'o'
            plt.plot(self.f, self.data.Tb[j], label=str(b))

        plt.plot(self.data.f, Tb, lt)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Brightness temperature [K]')


# ##############################################################################################################
#                                           ATMOSPHERE
# ##############################################################################################################
class atm_plots:
    def __init__(self, atm):
        self.atm = atm

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

    def plotTP(self, plot='auto'):
        """Plot the T-P profile"""
        if plot == 'auto':
            plt.figure(self.atm.planet + ': T-P')
        plt.title(self.atm.planet + ':  T-P profile')
        plt.loglog(self.atm.gas[self.atm.config.C['T']], self.atm.gas[self.atm.config.C['P']], label='T')
        self._frame_plot('T [K]', show_legend=False)

    def plotCloud(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
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

    def plotGas(self, dontPlot=['Z', 'P', 'T', 'DZ'], plot='auto'):
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

    def plotProp(self, dontPlot=['Z', 'P', 'T'], plot='auto'):
        if plot == 'auto':
            plt.figure(self.atm.planet + ': Properties')
        plt.title(self.atm.planet + ': other')
        for other in self.atm.config.LP:
            present, g = self.atm.is_present(self.atm.layerProperty[self.atm.config.LP[other]])
            if other in dontPlot or not present:
                continue
            plt.loglog(g, self.atm.gas[self.atm.config.C['P']], label=other)
        self._frame_plot('Property Value')


# ##############################################################################################################
#                                           SHAPE
# ##############################################################################################################
class shape_plots:
    def __init__(self, shp):
        self.shp = shp

    def plotShapes(self, planet, r, lat=90.0, delta_lng=0.0, gtypes=['ellipse', 'sphere', 'reference', 'gravity'], latstep='default'):
        colors = ['k', 'r', 'g', 'b', 'y', 'm', 'c']
        plt.figure('Shapes')
        for i, gtype in enumerate(gtypes):
            self.saveShape = [['oui']]
            rmag = self.shp.calcShape(planet, r, pclat=90.0, delta_lng=delta_lng, gtype=gtype, latstep=latstep)
            del self.shp.saveShape[0]  # remove nonsense first term
            self.shp.saveShape = np.array(self.shp.saveShape)
            if gtype != 'gravity':
                self.shp.saveShape = np.flipud(self.shp.saveShape)
            _y = []
            _z = []
            for v in self.shp.saveShape:
                _y.append(v[_Y])
                _z.append(v[_Z])
            plt.plot(_z, _y, color=colors[self.colorCounter % len(colors)], label=gtype)
            self.colorCounter += 1
            plt.axis('image')
            plt.legend()
