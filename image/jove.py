import pyfits
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import math
import cmath
import datetime
import time
import os
import string
import os


import img as __img__
import utils


# import aipy
stokes = {'rr': 0, 'll': 1}
units = {'natural': 1.0,
         'seconds': 1.0, 'minutes': 60.0, 'hours': 3600.0, 'days': 86400.0,
         'Hz': 1.0, 'MHz': 1.0E6, 'GHz': 1.0E9,
         'm': 1.0, 'cm': 0.1,
         'radians': 1.0, 'mas': 4.848136812E-9, 'degrees': 1.745329252E-2}
flagTypeKeys = {'reset': 0, 'time': 1, 'antenna': 2, 'baseline': 3, 'visibility': 4, 'phaseClosure': 5, 'weight': 6}
headerLineIndicator = '#$!^&'
omegaE = 7.29115E-5
# ##Generic Functions


def listFiles(ftype='fits', directory='.'):
    Files = []
    for line in os.listdir(directory):
        if ftype == 'all':
            Files.append(line)
        else:
            fil = line.split('.')
            if len(fil) > 1:
                if ftype == 'ripl':
                    ftypes = ['fits', 'image']  # ,'map','beam','lm']
                else:
                    ftypes = [ftype]
                for ft in ftypes:
                    if fil[-1].lower() == ft:
                        Files.append(line)
                        break
    for i, f in enumerate(Files):
        print '\t' + str(i) + '  ' + f
    return Files


ls = listFiles('ripl')


def convert2DMSfS(v):
    """Converts v to DMS"""
    D = int(v)
    v1 = 60.0 * (v - D)
    M = int(v1)
    v2 = 60.0 * (v1 - M)
    S = int(v2)
    fS = v2 - S
    return D, M, S, fS


def printPrettyDMS(v):
    dms = convert2DMSfS(v)
    pp = string.zfill(str(dms[0]), 2) + '  ' + string.zfill(str(dms[1]), 2) + '  ' + string.zfill(str(dms[2]), 2) + str(dms[3])[1:8]
    return pp


def julian2datetime(jultime, MJD=False):
    if not MJD:
        jultime -= 2400000.5
    jul0 = datetime.datetime(1858, 11, 17, 0, 0, 0)
    d = int(jultime)
    h1 = 24.0 * (jultime - d)
    h = int(h1)
    m1 = 60.0 * (h1 - h)
    m = int(m1)
    s1 = 60.0 * (m1 - m)
    s = int(s1)
    juldt = datetime.timedelta(days=d, hours=h, minutes=m, seconds=s)
    jul = jul0 + juldt
    return jul


def datetime2julian(dati, MJD=False):
    jul0 = datetime.datetime(1858, 11, 17, 0, 0, 0)
    jul = datetime.timedelta.total_seconds(dati - jul0) / units['days']
    if not MJD:
        jul += 2400000.5
    return jul


class Img():
    """Reads FITS files to process planet image data"""

    def __init__(self, srcFile='NEP-X-FINAL_KELVIN.FITS'):
        """Initialization opens the fits file and reads it in"""
        self.srcFile = srcFile

        # ##Some random initializations
        self.radecCtrArchive = [[], []]
        self.applyFlags = False
        self.flagDescriptor = []
        print datetime.datetime.now()

        # ## Open file and fix data column names (which doesn't work for some reason)
        self.p = pyfits.open(srcFile)
        self.p.info()
        # ##Get desired header info
        objName = self.p[0].header.get('OBJECT')
        print '\nObject:  ' + objName
        freqGHz = self.p[0].header.get('CRVAL3') / 1.0E9
        band = utils.getRFband(freqGHz, 'GHz')
        print 'Frequency:  %.4f  (%s band)' % (freqGHz, band)
        dateObs = self.p[0].header.get('DATE-OBS')
        dateMap = self.p[0].header.get('DATE-MAP')
        print 'Observed on ' + dateObs
        print 'Map made ' + dateMap
        self.res = [[], []]
        self.bw = [[], []]
        for i in range(len(self.p[0].header)):
            if type(self.p[0].header[i]) == str:
                hdata = self.p[0].header[i].split()
                if len(hdata) > 0:
                    if hdata[0] == 'IMAGR':
                        if hdata[1] == 'CELLSIZE(':
                            cellNo = int(hdata[2].strip(')')) - 1
                            v = self.p[0].header[i].split('=')
                            self.res[cellNo] = float(v[1])
                        if hdata[1] == 'BMAJ':
                            v = self.p[0].header[i].split('=')
                            self.bw[0] = float(v[1])
                        if hdata[1] == 'BMIN':
                            v = self.p[0].header[i].split('=')
                            self.bw[1] = float(v[1])
        print 'Resolution:  %.2f x %.2f' % (self.res[0], self.res[1])
        print 'Beamwidth:  %.2f x %.2f' % (self.bw[0], self.bw[1])

        # ##Set up table and field dictionaries
        self.fTab = {}
        for i in range(len(self.p)):
            nm = self.p[i].name.strip()
            if len(nm) == 0:
                nm = 'data'
            self.fTab[nm] = i
        # self.data = np.flipud(self.p[0].data[0][0])  # flip to convert fits orientation to standard image
        self.data = self.p[0].data[0][0]
        self.xyextents = [-self.res[0] * len(self.data) / 2.0, self.res[0] * len(self.data) / 2.0,
                          -self.res[1] * len(self.data) / 2.0, self.res[1] * len(self.data) / 2.0]
        self.x = []
        self.y = []
        for i in range(len(self.data)):
            self.x.append(self.xyextents[0] + i * self.res[0])
            self.y.append(self.xyextents[2] + i * self.res[0])
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        self.show()

    def show(self):
        __img__.showImage(self.x, self.y, self.data)

    def shift(self, xp, yp):
        """Just shift x, y values by (xp,yp) pixels to shift image location"""
        self.x = []
        self.y = []
        for i in range(len(self.data)):
            self.x.append(self.xyextents[0] + (i - xp) * self.res[0])
            self.y.append(self.xyextents[2] + (i - yp) * self.res[0])
        self.x = np.array(self.x)
        self.y = np.array(self.y)


class uv():
    """Reads FITS files to process planet image data"""

    def __init__(self, srcFile='NEP-X-FINAL_KELVIN.FITS'):
        """Initialization opens the fits file and reads it in"""
        self.srcFile = srcFile
        print datetime.datetime.now()

        # ##Open file and fix data column names (which doesn't work for some reason)
        self.p = pyfits.open(srcFile)
        self.p.info()

        # ##Set up table and field dictionaries
        self.fTab = {}
        self.hdr = []
        for i in range(len(self.p)):
            self.hdr.append(self.p[i].header)
            nm = self.p[i].name.strip()
            if len(nm) == 0:
                nm = 'data'
            self.fTab[nm] = i

        uvtab = self.fTab['AIPS UV']
        objName = self.p[uvtab].header.get('OBJECT')
        print '\nObject:  ' + objName
        freqGHz = self.p[uvtab].header.get('3CRVL6') / 1.0E9
        print 'Frequency:  ', freqGHz
        dateObs = self.p[uvtab].header.get('DATE-OBS')
        dateMap = self.p[uvtab].header.get('DATE-MAP')
        print 'Observed on ' + dateObs
        self.uu = []
        self.vv = []
        self.ww = []
        self.day = []
        numVis = len(self.p[uvtab].data)
        print 'Reading %d visibilities' % numVis
        for i in range(numVis):
            self.uu.append(self.p[uvtab].data[i][0])
            self.vv.append(self.p[uvtab].data[i][1])
            self.ww.append(self.p[uvtab].data[i][2])
            self.day.append(self.p[uvtab].data[i][4])

        totTime = 24.0 * (self.day[-1] - self.day[0])
        print 'Total time of observation is %f hours' % totTime
