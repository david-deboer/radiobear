import matplotlib.pyplot as plt
import numpy as np
import os.path
from radiobear import fileIO
from radiobear import utils


# ##############################################################################################################
#                                          GENERAL FILE PLOTTING
# ##############################################################################################################
def Tb(fn=None, xaxis='Frequency', directory='Output', file_type='spectrum', **kwargs):
    """plots brightness temperature against frequency and disc location:
           fn = filename to read (None will search...)
           xaxis = 'f[requency]' | 'w[avelength' ['freq']
           directory = subdirectory for data
           kwargs options are:
            xlog:  plot x-axis on log-scale if True
            ylog:  plot y-axis on log-scale if True
            legend:  include legend if True
           """

    fio = fileIO.FileIO(directory=directory)
    fio.read(fn=fn, file_type=file_type)

    # Frequency plot
    plt.figure('TB')
    for filen in fio.files:
        for i, b in enumerate(fio.data[filen].b):
            if xaxis[0].lower() == 'f':
                plotx = fio.data[filen].f
                xlabel = 'Frequency [GHz]'
            else:
                plotx = (utils.speedOfLight / 1E7) / fio.data[filen].f
                xlabel = 'Wavelength [cm]'
            label = "{}: {}".format(b, os.path.basename(filen))
            plt.plot(plotx, fio.data[filen].Tb[i], label=label)
            if 'xlog' in kwargs.keys() and kwargs['xlog']:
                plt.xscale('log')
            if 'ylog' in kwargs.keys() and kwargs['ylog']:
                plt.yscale('log')
    plt.xlabel(xlabel)
    plt.ylabel('Brightness Temperature [K]')
    if 'legend' in kwargs.keys() and kwargs['legend']:
        plt.legend()
    return fio


def Obs(fn, cols=[0, 1, 2], color='b', marker='o', delimiter=None, comline='!'):
    """
    <<<<????>>>>
    """
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
