import matplotlib.pyplot as plt
import numpy as np


# ##############################################################################################################
#                                          GENERAL FILE PLOTTING
# ##############################################################################################################
def TB(fn=None, xaxis='Frequency', xlog=False, justFreq=False, directory='Output', distance=4377233696.68):
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


def Obs(fn, cols=[0, 1, 2], color='b', marker='o', delimiter=None, comline='!'):
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
