import matplotlib.pyplot as plt
import numpy as np


# ##############################################################################################################
#                                              BRIGHTNESS
# ##############################################################################################################
class plots:
    def __init__(self, bright):
        self.bright = bright

    def show(self):
        plt.show()

    def alpha(self):
        # ####-----Alpha
        plt.figure('alpha')
        for i, f in enumerate(self.bright.freqs):
            if self.bright.output_type == 'frequency':
                label = ('a: {:.1f} GHz').format(f)
            else:
                label = ('{:.1f} cm').format(30.0 / f)
            plt.loglog(self.bright.layerAlpha[i][1:], self.bright.P, label=label)
        self._frame_plot('dB/km', True)

    def intW(self):
        # ####-----Weigthing functions
        plt.figure('INT_W')
        lt = '-'
        if len(self.bright.freqs) == 1:
            lt = 'o'
        plt.plot(self.bright.freqs, self.bright.integrated_W, lt)
        plt.title('Integrated weighting function')
        plt.xlabel('Frequency [GHz]')

    def W(self, normW4plot):
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

    def raypath(self):
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

    def observer(self, b, req=None, rpol=None):
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
