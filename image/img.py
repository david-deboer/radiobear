from __future__ import absolute_import, division, print_function
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import numpy as np
import math
from scipy import misc
from scipy import interpolate
import prog_path
import utils
import shape
import raypath
import fileIO


class Img:
    def __init__(self, i='?', readAndPlot=True):
        self.scalar = (float, int)
        self.TB = fileIO.FileIO('frequency')
        self.kerneled = False
        self.kernel = None
        self.convolved = False
        self.reducedImagesMade = False
        if readAndPlot:
            self.getImage(i)

    def getImage(self, ifilename=None, usedir='Output'):
        """Reads in image files from pyPlanet"""

        self.TB.read(ifilename)  # what to do about usedir?
        print('For now copy data over to local self...')
        self.data = self.TB.data
        self.header = self.TB.header
        self.xyextents = self.TB.xyextents
        self.x = self.TB.x
        self.y = self.TB.y
        self.resolution = self.TB.resolution
        self.show()

    def fix(self, fixThreshhold='auto', padValue='auto'):
        """Fix or reset points:
              padValue = float or 'auto'=nearest neighbor average
              fixThreshhold = float sets all points above trim to padValue
              fixThreshhold = list of pairs, sets those points to padValue
              fixThreshhold = 'auto', set to 110% of non-zero average"""
        print('called with: ', fixThreshhold)
        fixVal = 0
        if type(fixThreshhold) == list:
            if len(fixThreshhold) == 2:
                fixVal = fixThreshhold
            else:
                print('Error in ', fixThreshhold)
                return None
        elif fixThreshhold[0].lower() == 'v':
            figname = fixThreshhold + str(padValue)
            plt.figure(figname)
            for i in range(len(self.data)):
                plt.plot(self.data[i])
            return 0
        elif fixThreshhold[0].lower() == 'd':
            print("Derivative fix...")
            ft = self.fixDeriv(padValue)
            self.fix(ft)
            return '...done'
        elif fixThreshhold == 'auto':
            dataShape = np.shape(self.data)
            nave = 0
            ave = 0.0
            for i in range(dataShape[0]):
                for j in range(dataShape[1]):
                    if self.data[i, j] > 0.0:
                        nave += 1
                        ave += self.data[i, j]
            ave = ave / nave
            fixThreshhold = 1.1 * ave
            print('Non-zero average = {:.2f}, set fixThreshhold = {:.2f}'.format(ave, fixThreshhold))
            fixVal = np.where(self.data > fixThreshhold)
        else:
            fixVal = np.where(self.data > fixThreshhold)

        if padValue == 'auto':
            offset = [[1, -1], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1]]
            padValue = []
            for i in range(len(fixVal[0])):
                padVal = 0.0
                for oval in offset:
                    padVal += self.data[fixVal[0][i] + oval[0], fixVal[1][i] + oval[1]]
                padValue.append(padVal / len(offset))
        else:
            padVal = padValue
            padValue = []
            for i in range(len(fixVal[0])):
                padValue.append(padVal)

        print('Pad values')
        for i in range(len(fixVal[0])):
            print('point {}, {} = {:.2f}'.format(fixVal[0][i], fixVal[1][i], padValue[i]))
            self.data[fixVal[0][i], fixVal[1][i]] = padValue[i]
        return len(fixVal[0])

    def fixDeriv(self, threshhold=140.0):
        """Fixes points based on derivative in planet disc"""
        if threshhold == 'auto':
            print('You probably forgot to provide a value for fixDeriv in fix...')
            return 0
        numFix = 0
        fixLoc = [[], []]
        plt.figure('fixDeriv')
        for i in range(len(self.data)):
            inDisc = False
            for j in range(len(self.data[i])):
                if self.data[i, j] > 0.0 and not inDisc:
                    inDisc = True
                elif inDisc and self.data[i, j] < 1E-6:
                    inDisc = False
                if not inDisc:
                    continue

                derivForward = self.data[i, j + 1] - self.data[i, j]
                derivReverse = self.data[i, j] - self.data[i, j - 1]
                if abs(derivForward) > threshhold:
                    print(i, j, derivForward)
                    plt.plot(self.data[i])
                    plt.plot(j, self.data[i, j], 'o')
                    numFix += 1
                    fixLoc[0].append(i)
                    fixLoc[1].append(j)
                elif abs(derivReverse) > threshhold:
                    print(i, j, derivReverse)
                    plt.plot(self.data[i])
                    plt.plot(j - 1, self.data[i, j], 'o')
                    numFix += 1
                    fixLoc[0].append(i)
                    fixLoc[1].append(j)

        print('Number found:  {}'.format(len(fixLoc[0])))
        return fixLoc

    def outline(self, tip='auto', reference='ref'):
        """Finds approximate elliptical outline.
              - tip:  angle in the plane of the sky, 'auto' to match observation, 0 for aligned vertical
              - reference:  ['edge' or 'ref' or number_km]"""
        outline = [[], []]
        try:  # this is all done in km
            radii = []
            for i in range(len(self.header['radii'])):
                if isinstance(self.header['radii'][i], self.scalar):
                    radii.append(self.header['radii'][i])
            Req = np.max(radii)
            Rpol = np.min(radii)
            rNorm = self.header['rNorm'][0]
            dist = self.header['distance'][0] * utils.Units[self.header['distance'][1]] / utils.Units['km']
            checkResolution = (180.0 / math.pi) * 3600.0 * math.atan(self.header['b'][0] * rNorm / dist)
            lat_pg = self.header['orientation'][1] * np.pi / 180.0
        except KeyError:
            print('Not enough header information')
            return None

        print('Reference disc is {:.1f} x {:.1f} km at {:.1f} km distance'.format(Req, Rpol, dist))
        print('Image resolution is %f arcsec' % (self.resolution))
        if isinstance(reference, self.scalar):
            a = float(reference)
        elif reference == 'edge':
            a = self.header['rNorm'][0]
        else:  # assume reference (1-bar) level
            a = Req
        # print 'assuming reference diameter of %.1f' % (a)
        f = (Req - Rpol) / Req
        lat_pc = math.atan(math.tan(lat_pg) * (1.0 - f)**2)
        b = a * (1.0 - f * math.cos(lat_pc)**2)

        # ## Convert to arcsec
        a = (180.0 / math.pi) * 3600.0 * math.atan(a / dist)
        b = (180.0 / math.pi) * 3600.0 * math.atan(b / dist)
        for phi in np.arange(0.0, 2.0 * np.pi, np.pi / 180.0):
            outline[0].append(a * math.cos(phi))
            outline[1].append(b * math.sin(phi))
        if tip == 'auto':
            tip = self.header['orientation'][0]
        if tip:
            tip = tip * np.pi / 180.0
            for i in range(len(outline[0])):
                r = np.array([outline[0][i], outline[1][i], 0.0])
                r = shape.rotZ(tip, r)
                outline[0][i] = r[0]
                outline[1][i] = r[1]
        return outline

    def lat(self, lat_pg, tip='auto', rotate='auto', lngs=360.0, reference='ref', verbose=False):
        """Calculates the projected latitude on the disc.
               - lat_pg:  desired planetographic latitude in degrees.
               - tip:  used to tip to orientation ('auto') or another specific value (e.g. 0.0 to keep vertical)
               - rotate:  used to rotate for sub-earth point 'auto' matches the header value
               - reference:  ['edge' or 'ref' (usually 1-bar) or number_km]
               - lngs:  angle around for the latitude line - somewhat trial and error."""
        try:   # ## in km
            radii = []
            for i in range(len(self.header['radii'])):
                if isinstance(self.header['radii'][i], self.scalar):
                    radii.append(self.header['radii'][i])
            Req = np.max(radii)
            Rpol = np.min(radii)
            rNorm = self.header['rNorm'][0]
            dist = self.header['distance'][0] * utils.Units[self.header['distance'][1]] / utils.Units['km']
            checkResolution = (180.0 / math.pi) * 3600.0 * math.atan(self.header['b'][0] * rNorm / dist)
        except KeyError:
            print('Not enough header information')
            return None
        f = (Req - Rpol) / Req
        lat_pg = lat_pg * math.pi / 180.0
        lat_pc = math.atan(math.tan(lat_pg) * (1.0 - f)**2)
        if isinstance(reference, self.scalar):
            a = float(reference)
        elif reference == 'edge':
            a = self.header['rNorm'][0]
        else:  # assume reference (1-bar) level
            a = Req
        # print 'assuming reference diameter of %.1f' % (a)
        b = (1.0 - f) * a

        # ## Convert to arcsec
        a = (180.0 / math.pi) * 3600.0 * math.atan(a / dist)
        b = (180.0 / math.pi) * 3600.0 * math.atan(b / dist)

        ylat = b * math.sin(lat_pc)
        alat = a * math.cos(lat_pc)
        if tip == 'auto':
            tip = self.header['orientation'][0]
        tip *= math.pi / 180.0
        if rotate == 'auto':
            rotate = self.header['orientation'][1]
        rotate = math.atan(math.tan(rotate * math.pi / 180.0) * (1.0 - f)**2)
        if verbose:
            print('f, a, b = ', f, a, b)
            print('lat_pg --> lat_pc = ', 180.0 * lat_pg / np.pi, 180.0 * lat_pc / np.pi)
            print('ylat, alat = ', ylat, alat)
            print('tip, rotate = ', 180.0 * tip / np.pi, 180.0 * rotate / np.pi)
        lat = [[], []]
        for th in np.arange(-lngs / 2.0, lngs / 2.0, lngs / 200.0):
            angle = th * np.pi / 180.0
            r = np.array([alat * math.sin(angle), ylat, alat * math.cos(angle)])
            r = raypath.__rotate2obs__(rotate, tip, r)
            lat[0].append(r[0])
            lat[1].append(r[1])
        return lat

    def drawLat(self, tip='auto', reference='ref'):
        """ad hoc function to draw specific latitude lines"""
        out = self.outline(reference=reference)
        plt.plot(out[0], out[1], 'k')
        lat = self.lat(-89.0, rotate='auto', tip=tip, reference=reference)
        plt.plot(lat[0], lat[1], 'k')
        lat = self.lat(-60.0, lngs=260.0, rotate='auto', tip=tip, reference=reference)
        plt.plot(lat[0], lat[1], 'k--')
        lat = self.lat(-30.0, lngs=200.0, rotate='auto', tip=tip, reference=reference)
        plt.plot(lat[0], lat[1], 'k--')
        lat = self.lat(0.0, lngs=192.0, rotate='auto', tip=tip, reference=reference)
        plt.plot(lat[0], lat[1], 'k')
        lat = self.lat(30.0, lngs=160.0, rotate='auto', tip=tip, reference=reference)
        plt.plot(lat[0], lat[1], 'k--')

    def saveImage(self, header=None, image=None, filename=None):
        """Saves header/image to file"""
        if image is None:
            image = self.data
            header = self.header
        if filename is None:
            print('Saving image composed of file headers:')
            print(self.header)
            filename = raw_input('output filename: ')
        fp = open(filename, 'w')
        for hkey in header.keys():
            s = '#  ' + hkey + ':  '
            for v in header[hkey]:
                s += str(v) + '  '
            fp.write(s + '\n')
            print(s)
        for row in image:
            s = ''
            for col in row:
                s += '%.2f\t' % (col)
            fp.write(s[0: -1] + '\n')
        fp.close()

    def editImageB(self, bmod, Tmod):
        import planet
        nedit = 0
        p = planet.planet('functions')
        bdata = p.bRequest(self.header['b'][0],block=[1,1])
        if len(bdata) != np.size(self.data):
            print('b-values (%d) do not number the same as the image size (%d)' % (len(bdata),np.size(self.data)))
            return 0
        if len(bmod)==3:    # assume center and radius
            print('Circle edit')
            pixel_modlist = []
            for i in range(len(bdata)):
                if math.sqrt( bdata[i][0]**2 + bdata[i][1]**2 ) > 0.95:
                    continue  # this is off the disk (assumed circular...)
                elif math.sqrt( (bdata[i][0] - bmod[0])**2 + (bdata[i][1] - bmod[1])**2) < bmod[2]:
                    row = int(float(i) // float(len(self.data)))
                    col = i - row*len(self.data)
                    pixel_modlist.append([row,col])
                    nedit+=1
        elif len(bmod)==4:  # assume square
            print('Square')
        else:               # assume list
            print('list')
            if type(Tmod) == float:
                tmp = Tmod
                Tmod = []
                for b in bmod:
                    Tmod.append(tmp)
        for px in pixel_modlist:
            self.data[px[0],px[1]] = Tmod
        return nedit

    def editImageLat(self,latmod,Tmod,limb=[0.0]):
        """This only works for negative latitudes to pole.  Need to make northern,southern,straddle + bands.
            latmod is the planetographic latitude that marks the boundary (degrees)
            Tmod is/are the temperature values
            limb are the corresponding b-values"""
        import planet
        import atmosphere
        pyPlanetPath = os.getenv('PYPLANETPATH')
        if pyPlanetPath == None:
            print('No PYPLANETPATH environment variable')
            pyPlanetPath = './'
        atm = atmosphere.atmosphere('neptune',pyPlanetPath)
        atm.run()
        p = planet.planet('functions')
        bdata = p.bRequest(self.header['b'][0],block=[1,1])
        if len(bdata) != np.size(self.data):
            print('b-values (%d) do not number the same as the image size (%d)' % (len(bdata),np.size(self.data)))
            return 0
        radii = []
        for i in range(len(self.header['radii'])):
            if isinstance(self.header['radii'][i], self.scalar):
                radii.append(self.header['radii'][i])
        rNorm = self.header['rNorm'][0]
        Req = np.max(radii)
        Rpol = np.min(radii)
        f = (Req-Rpol)/Req
        gtype = 'ellipse'
        tip = self.header['aspect'][0]
        tip*=math.pi/180.0
        rotate = self.header['aspect'][1]
        #rotate = -1.0*math.atan(math.tan(rotate*math.pi/180.0)*(1.0-f)**2)
        rotate*=math.pi/180.0
        lat_pg = latmod*math.pi/180.0
        lat_pcmod = math.atan(math.tan(lat_pg)*(1.0-f)**2)
        pixel_modlist = []
        nedit = 0
        if isinstance(Tmod, self.scalar):
            Tmod = list(Tmod)
        if len(Tmod) != len(limb):
            print('limb is incorrectly specified')
            return 0
        for i in range(len(bdata)):
            edge, bmod = raypath.__findEdge__(atm,bdata[i],rNorm,tip,rotate,gtype,printdot=False)
            if edge == None:
                continue
            pclat = math.asin( np.dot(edge,raypath.yHat)/np.linalg.norm(edge) )
            delta_lng = math.atan2( np.dot(edge,raypath.xHat), np.dot(edge,raypath.zHat) )
            if pclat > 0.0:
                break
            if pclat < lat_pcmod:
                row = int(float(i) // float(len(self.data)))
                col = i - row*len(self.data)
                pixel_modlist.append([row,col,i])
                nedit+=1
        breakb = 0.96
        if len(Tmod) > 1:
            wherebb = np.where(np.array(limb) < breakb)
            Tfunc3 = interpolate.interp1d(limb[wherebb[0][0]:wherebb[0][-1]+1],Tmod[wherebb[0][0]:wherebb[0][-1]+1],kind=3,bounds_error=False,fill_value=0.0)
            Tfunc1 = interpolate.interp1d(limb,Tmod,kind=1)
        else:
            Tfunc3 = lambda x: Tmod[0]
            Tfunc1 = lambda x: Tmod[0]
        testPlotLimb = False
        if testPlotLimb:
            pltl = np.linspace(limb[0],limb[-1],50)
            pltt = Tfunc3(pltl)
            plt.figure('limbFunction')
            plt.plot(pltl,pltt)
            pltt = Tfunc1(pltl)
            plt.plot(pltl,pltt)
            plt.plot(limb,Tmod,'o')
            math.sqrt(-1.)
        for px in pixel_modlist:
            bval = math.sqrt( bdata[px[2]][0]**2 + bdata[px[2]][1]**2 )
            if bval < limb[0]:
                bval = limb[0]
            elif bval > limb[-1]:
                bval = limb[-1]
            if bval < breakb:
                Tnew = Tfunc3(bval)
            else:
                Tnew = Tfunc1(bval)
            self.data[px[0],px[1]] = Tnew
        return nedit

    def generateKernel(self, fwhm, ktype='gaussian'):
        """Generate the convolving kernel.  ktype=gaussian is only option right now."""
        self.kernel = np.zeros(np.shape(self.data))
        alpha = 4.0*math.log(2)/fwhm**2.0
        self.fwhm = fwhm
        res = self.resolution
        print('Generating '+ktype+' kernel with fwhm = ',fwhm)
        if type(res) == str:
            self.resolution = float(raw_input('input image resolution in arcsec:  '))
        center = [len(self.kernel)/2, len(self.kernel[0])/2]
        for i in range(len(self.kernel)):
            for j in range(len(self.kernel)):
                d = math.sqrt((float(i)-center[0])**2.0 + (float(j)-center[1])**2.0)*res
                self.kernel[j,i] = math.exp(-alpha*(d**2))
        self.header['kernel'] = [fwhm,'arcsec',ktype]
        print('kernel at self.kernel and self.kernelHeader')
        self.kerneled = True

        return alpha

    def discAverage(self,threshhold=1.0):
        nave = 0
        Tave = 0.0
        for row in self.data:
            for T in row:
                if T>threshhold:
                    nave+=1
                    Tave+=T
        Tave = Tave/nave
        print('Disc average = %f' % (Tave))
        return Tave


    def convolve(self, fwhm=None):
        """Convolves self.data with self.kernel (from self.generateKernel).  Output to self.convolvedImage"""
        if not self.kerneled and fwhm==None:
            print('Error:  no kernel')
            return None
        if fwhm:
            self.generateKernel(fwhm,ktype='gaussian')
        self.convolvedHeader = self.header
        self.header['note'] = ['convolved']
        self.convolvedImage = convolvend(self.data,self.kernel,normalize_kernel=True)
        print('convolved image at self.convolvedImage and self.convolvedHeader')
        self.convolved = True

    def rotate_image(self,image,rotate=None):
        """rotates and overwrites image.
               - rotate is angle in degrees"""
        if rotate == None:
            rotate = self.header['aspect'][0]
        print('Rotate:  ',rotate)
        img_max = np.max(image)
        if rotate != 0.0:
            image = misc.imrotate(image,rotate)
            maxout = float(np.max(image))
            image = img_max*image/maxout
        return image

    def diff(self, x, y, image, region='auto', mask='auto', resize='up', bestFit=True, kind='cubic', showMask=False, plot=False):
        """Diff'ing images:
              x - x data of observation ('n.x')
              y - y data of observation ('n.y')
              image - observation ('n.data')
              region - part of image to use (hard-coded right now!)
              mask = part of image to not use in fit (e.g. for the hot spot - hard-coded right now!)
              return  image - self.convolvedImage
              resize is 'up' or 'down'
              bestFit is True/False or a list.  len=2 computes at one offset, len=4 computes at those ranges"""

        if self.convolved==False:
            print('Image not convolved')
            return None
        if region == 'auto':  # for now
            if self.header['band'] == 'C':
                region = [-2.0,2.0,-2.0,2.0]
            elif self.header['band'] == 'Q':
                region = [-1.5,1.5,-1.5,1.5]
            else:
                region = [-1.5,1.5,-1.5,1.5]
        if mask == 'auto':    # for now
            mask = [-0.5,0.1,-0.95,-0.6]
        maskShown = False
        res1 = x[1] - x[0]
        res2 = self.resolution
        if not self.reducedImagesMade:
            self.reducedImagesMade = resize
            if (res1 < res2 and resize=='up') or (res2 < res1 and resize=='down'):
                print('Resizing '+resize+' to input image over region '+str(region))
                v2x = np.where( (self.x > region[0]) & (self.x < region[1]) )
                v2y = np.where( (self.y > region[2]) & (self.y < region[3]) )
                x2reduced = x[v2x]
                y2reduced = y[v2y]
                image2reduced = self.convolvedImage[v2x[0][0]:v2x[0][-1]+1,v2y[0][0]:v2y[0][-1]+1]
                cbef = str(np.shape(image2reduced))
                resFunc = interpolate.interp2d(x2reduced,y2reduced,image2reduced,kind=kind,fill_value=0.0)
                v1x = np.where( (x > region[0]) & (x < region[1]) )
                v1y = np.where( (y > region[2]) & (y < region[3]) )
                x1reduced = x[v1x[0]]
                y1reduced = y[v1y[0]]
                self.image2reduced = self.convolvedImage[v1x[0][0]:v1x[0][-1],v1y[0][0]:v1y[0][-1]]
                caft = str(np.shape(self.image2reduced))
                self.image1reduced = resFunc(x1reduced,y1reduced)
                ibef = iaft = str(np.shape(self.image2reduced))
                self.xdiff = x1reduced
                self.ydiff = y1reduced
                del image2reduced
            elif (res1 < res2 and resize=='down') or (res2 < res1 and resize=='up'):
                print('Resizing '+resize+' to convolvedImage over region '+str(region))
                v1x = np.where( (x > region[0]) & (x < region[1]) )
                v1y = np.where( (y > region[2]) & (y < region[3]) )
                x1reduced = x[v1x]
                y1reduced = y[v1y]
                image1reduced = image[v1x[0][0]:v1x[0][-1]+1,v1y[0][0]:v1y[0][-1]+1]
                ibef = str(np.shape(image1reduced))
                resFunc = interpolate.interp2d(x1reduced,y1reduced,image1reduced,kind=kind,fill_value=0.0)
                v2x = np.where( (self.x > region[0]) & (self.x < region[1]) )
                v2y = np.where( (self.y > region[2]) & (self.y < region[3]) )
                x2reduced = self.x[v2x[0]]
                y2reduced = self.y[v2y[0]]
                self.image2reduced = self.convolvedImage[v2x[0][0]:v2x[0][-1]+1,v2y[0][0]:v2y[0][-1]+1]
                self.image1reduced = resFunc(x2reduced,y2reduced)
                iaft = str(np.shape(self.image1reduced))
                cbef = caft = str(np.shape(self.image2reduced))
                self.xdiff = x2reduced
                self.ydiff = y2reduced
                del image1reduced
            else:
                print('screwed up somewhere')
                return None
            print('\tshape(image):  %s -> %s' % (ibef,iaft))
            print('\tshape(convolved):  %s -> %s' % (cbef,caft))
            if plot:
                plt.figure(101)
                plt.subplot(221)
                showImage(self.xdiff,self.ydiff,self.image1reduced)
                plt.subplot(222)
                showImage(self.xdiff,self.ydiff,self.image2reduced)
        if showMask:
            maskImg = np.copy(self.image1reduced)

        if bestFit:    # check for best fit, assume convolvedImage is the centered one
            if type(bestFit) == list:
                if len(bestFit) == 2:
                    col_offset_list = [bestFit[0]]
                    row_offset_list = [bestFit[1]]
                elif len(bestFit) == 4:
                    col_offset_list = range(bestFit[0],bestFit[1]+1)
                    row_offset_list = range(bestFit[2],bestFit[3]+1)
                else:
                    print('Error in bestFit list')
                    col_offset_list = [0]
                    row_offset_list = [0]
            else:
                col_offset_list = range(-10,11)
                row_offset_list = range(-10,11)
            print('Finding best residual over range ',)
            print('('+str(col_offset_list[0])+','+str(row_offset_list[0])+') - ('+str(col_offset_list[-1])+','+str(row_offset_list[-1])+')')
            print('with mask '+str(mask))
            vMasky = np.where( (self.xdiff > mask[0]) & (self.xdiff < mask[1]) )
            vMaskx = np.where( (self.ydiff > mask[2]) & (self.ydiff < mask[3]) )
            if showMask and not maskShown:
                maskImg[vMaskx[0][0]:vMaskx[0][-1]+1,vMasky[0][0]:vMasky[0][-1]+1] = 0.0
                plt.figure('Mask')
                showImage(self.xdiff,self.ydiff,maskImg)
                maskShown = True
            best = 1.0E12
            col_best = 0
            row_best = 0
            for row_offset in row_offset_list:
                print('Trial:  ',row_offset)
                for col_offset in col_offset_list:
                    sqerr = 0.0
                    for i in range( len(self.image1reduced) ):
                        for j in range( len(self.image1reduced[0]) ):
                            if (i in vMaskx[0]) and (j in vMasky[0]):
                                sqerr+=0.0
                            else:
                                try:
                                    sqerr+= (self.image1reduced[i+row_offset,j+col_offset] - self.image2reduced[i,j])**2
                                except IndexError:
                                    sqerr+=0.0
                    sqerr/=np.size(self.image1reduced)
                    if sqerr < best:
                        best = sqerr
                        row_best = row_offset
                        col_best = col_offset
            best = math.sqrt(best)
            print('Offset:  %d, %d with residual %f' % (col_best,row_best, best))
            self.residual = best
        else:
            row_best = 0
            col_best = 0

        # difference images
        self.diffImage = np.zeros(np.shape(self.image1reduced))
        for i in range( len(self.image1reduced) ):
            for j in range( len(self.image1reduced[0]) ):
                try:
                    self.diffImage[i,j] = self.image1reduced[i+row_best,j+col_best] - self.image2reduced[i,j]
                except IndexError:
                    self.diffImage[i,j] = 0.0
        return self.diffImage

    def beam(self,pos=[-1.2,-1.2]):
        bm = [[],[]]
        for th in np.arange(-2.0*np.pi,2.0*np.pi,np.pi/180.0):
            bm[0].append(self.fwhm*math.cos(th)/2.0 + pos[0])
            bm[1].append(self.fwhm*math.sin(th)/2.0 + pos[1])
        plt.plot(bm[0],bm[1],'k')


    def show(self,imtype='r'):
        if imtype[0] == 'c':
            showImage(self.x,self.y,self.convolvedImage)
        elif imtype[0] == 'k':
            showImage(self.x,self.y,self.kernel)
        elif imtype[0] == 'd':
            showImage(self.xdiff,self.ydiff,self.diffImage)
        else:
            showImage(self.x,self.y,self.data)

    def diffmatrix(self,cols,rows,x,y,image):
        """This does a bunch of fits under diff and plots in an image-matrix format"""
        plt.figure('matrix')
        npc = len(rows)*len(cols)
        pc = 1
        for j in rows:
            for i in cols:
                diff = self.diff(x,y,image,bestFit=[i,j])
                plt.subplot(len(rows),len(cols),pc)
                print('%d of %d' % (pc,npc))
                pc+=1
                self.show(imtype='d')
                self.drawLat()
                plt.title('(%d, %d)  %.2f K' % (i,j,self.residual))



#################END OF CLASS################
def showImage(x,y,image,norm='lin'):
    """This assumes that the 'image' is in data mode and imshows a flipud version"""
    extent = [x[0],x[-1],y[0],y[-1]]
    show = np.flipud(image)
    if norm=='log':
        plt.imshow(show,extent=extent,norm=LogNorm)
    else:
        plt.imshow(show,extent=extent)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('K',rotation='horizontal')

def diffOLD(image1,res1,image2,res2,resize='up',bestFit=True,plot=False):
    """Diff'ing centered images with size(image1) > size(image2):
          return  image1 - resize(image2)
          resize is 'up' or 'down'"""

    img_max =  [np.max(image1), np.max(image2)]
    if abs(res1 - res2)/(res1 + res2) < 0.001:
        print('resolutions nearly equal - no resize')
    else:
        factor = res1/res2
        if factor < 1.0 and resize=='up':
            factor = 1.0/factor
        elif factor > 1.0 and resize=='down':
            factor = 1.0/factor
        if (res1 < res2 and resize=='up') or (res2 < res1 and resize=='down'):
            image2 = misc.imresize(image2,factor,interp='bilinear')
            maxOut = float(np.max(image2))
            image2 = img_max[1]*(image2/maxOut)
        elif (res1 < res2 and resize=='down') or (res2 < res1 and resize=='up'):
            image1 = misc.imresize(image1,factor,interp='bilinear')
            maxOut = float(np.max(image1))
            image1 = img_max[0]*(image1/maxOut)
        else:
            print('screwed up somewhere')
            return None
    print('factor = ',factor)
    offset = [-2,-1,0,1,2]
    extra = [len(image2)%2, len(image2[0])%2]
    center1 = [len(image1)/2,len(image1[0])/2]
    center2 = [len(image2)/2,len(image2[0])/2]

    if bestFit:
    # check for best fit

        best = 1.0E12
        i_best = 0
        j_best = 0
        for i_offset in offset:
            for j_offset in offset:
                sqerr = 0.0
                for i in range( -center2[0], center2[0]+extra[0] ):
                    for j in range( -center2[1], center2[1]+extra[1] ):
                        sqerr+= (image1[center1[0]+i+i_offset,center1[1]+j+j_offset] - image2[center2[0]+i,center2[1]+j])**2
                sqerr/=np.size(image2)
                if sqerr < best:
                    best = sqerr
                    i_best = i_offset
                    j_best = j_offset
        print('Offset:  %d, %d with residual %f' % (i_best,j_best, math.sqrt(best)))
    else:
        i_best = 0
        j_best = 0
    # difference images

    diffImage = np.zeros(np.shape(image2))
    for i in range( -center2[0], center2[0]+extra[0] ):
        for j in range( -center2[1], center2[1]+extra[1] ):
            diffImage[i+center2[0],j+center2[1]] = image1[center1[0]+i+i_best,center1[1]+j+j_best] - image2[center2[0]+i,center2[1]+j]
    if plot:
        if type(plot)==int:
            plt.figure(plot)
        else:
            plt.figure(43)
        fs = len(diffImage)
        plt.subplot(131)
        plt.imshow(image1)
        plt.axis([center1[0]-fs/2,center1[0]+fs/2,center1[1]+fs/2,center1[1]-fs/2])
        plt.colorbar()
        plt.title('Observation')
        plt.subplot(132)
        plt.imshow(image2)
        plt.colorbar()
        plt.title('Convolved model')
        plt.subplot(133)
        plt.imshow(diffImage)
        plt.colorbar()
        plt.title('Difference')
    return diffImage

#################################################################################
# Copied from
#   http://code.google.com/p/agpy/source/browse/trunk/AG_fft_tools/convolve_nd.py
#################################################################################
import warnings

try:
    import fftw3
    has_fftw = True

    def fftwn(array, nthreads=1):
        array = array.astype('complex').copy()
        outarray = array.copy()
        fft_forward = fftw3.Plan(array, outarray, direction='forward',
                flags=['estimate'], nthreads=nthreads)
        fft_forward.execute()
        return outarray

    def ifftwn(array, nthreads=1):
        array = array.astype('complex').copy()
        outarray = array.copy()
        fft_backward = fftw3.Plan(array, outarray, direction='backward',
                flags=['estimate'], nthreads=nthreads)
        fft_backward.execute()
        return outarray / np.size(array)
except ImportError:
    fftn = np.fft.fftn
    ifftn = np.fft.ifftn
    has_fftw = False
# I performed some fft speed tests and found that scipy is slower than numpy
# http://code.google.com/p/agpy/source/browse/trunk/tests/test_ffts.py However,
# the speed varied on machines - YMMV.  If someone finds that scipy's fft is
# faster, we should add that as an option here... not sure how exactly


def convolvend(array, kernel, boundary='fill', fill_value=0,
        crop=True, return_fft=False, fftshift=True, fft_pad=True,
        psf_pad=False, interpolate_nan=False, quiet=False,
        ignore_edge_zeros=False, min_wt=0.0, normalize_kernel=False,
        use_numpy_fft=not has_fftw, nthreads=1):
    """
    Convolve an ndarray with an nd-kernel.  Returns a convolved image with shape =
    array.shape.  Assumes image & kernel are centered.

    Parameters
    ----------
    array: `numpy.ndarray`
          Array to be convolved with *kernel*
    kernel: `numpy.ndarray`
          Will be normalized if *normalize_kernel* is set.  Assumed to be
          centered (i.e., shifts may result if your kernel is asymmetric)

    Options
    -------
    boundary: str, optional
        A flag indicating how to handle boundaries:
            * 'fill' : set values outside the array boundary to fill_value
                       (default)
            * 'wrap' : periodic boundary
    interpolate_nan: bool
        attempts to re-weight assuming NAN values are meant to be ignored, not
        treated as zero.  If this is off, all NaN values will be treated as
        zero.
    ignore_edge_zeros: bool
        Ignore the zero-pad-created zeros.  This will effectively decrease
        the kernel area on the edges but will not re-normalize the kernel.
        This parameter may result in 'edge-brightening' effects if you're using
        a normalized kernel
    min_wt: float
        If ignoring NANs/zeros, force all grid points with a weight less than
        this value to NAN (the weight of a grid point with *no* ignored
        neighbors is 1.0).
        If `min_wt` == 0.0, then all zero-weight points will be set to zero
        instead of NAN (which they would be otherwise, because 1/0 = nan).
        See the examples below
    normalize_kernel: function or boolean
        if specified, function to divide kernel by to normalize it.  e.g.,
        normalize_kernel=np.sum means that kernel will be modified to be:
        kernel = kernel / np.sum(kernel).  If True, defaults to
        normalize_kernel = np.sum

    Advanced options
    ----------------
    fft_pad: bool
        Default on.  Zero-pad image to the nearest 2^n
    psf_pad: bool
        Default off.  Zero-pad image to be at least the sum of the image sizes
        (in order to avoid edge-wrapping when smoothing)
    crop: bool
        Default on.  Return an image of the size of the largest input image.
        If the images are asymmetric in opposite directions, will return the
        largest image in both directions.
        For example, if an input image has shape [100,3] but a kernel with shape
      [6,6] is used, the output will be [100,6].
    return_fft: bool
        Return the fft(image)*fft(kernel) instead of the convolution (which is
        ifft(fft(image)*fft(kernel))).  Useful for making PSDs.
    fftshift: bool
        If return_fft on, will shift & crop image to appropriate dimensions
    nthreads: int
        if fftw3 is installed, can specify the number of threads to allow FFTs
        to use.  Probably only helpful for large arrays
    use_numpy_fft: bool
        Force the code to use the numpy FFTs instead of FFTW even if FFTW is
        installed

    Returns
    -------
    default: `array` convolved with `kernel`
    if return_fft: fft(`array`) * fft(`kernel`)
      * if fftshift: Determines whether the fft will be shifted before
        returning
    if not(`crop`) : Returns the image, but with the fft-padded size
        instead of the input size

    Examples
    --------
    >>> convolvend([1,0,3],[1,1,1])
    array([ 1.,  4.,  3.])

    >>> convolvend([1,np.nan,3],[1,1,1],quiet=True)
    array([ 1.,  4.,  3.])

    >>> convolvend([1,0,3],[0,1,0])
    array([ 1.,  0.,  3.])

    >>> convolvend([1,2,3],[1])
    array([ 1.,  2.,  3.])

    >>> convolvend([1,np.nan,3],[0,1,0], interpolate_nan=True)
    array([ 1.,  0.,  3.])

    >>> convolvend([1,np.nan,3],[0,1,0], interpolate_nan=True, min_wt=1e-8)
    array([  1.,  nan,   3.])

    >>> convolvend([1,np.nan,3],[1,1,1], interpolate_nan=True)
    array([ 1.,  4.,  3.])

    >>> convolvend([1,np.nan,3],[1,1,1], interpolate_nan=True, normalize_kernel=True, ignore_edge_zeros=True)
    array([ 1.,  2.,  3.])

    """


    # Checking copied from convolve.py - however, since FFTs have real &
    # complex components, we change the types.  Only the real part will be
    # returned!
    # Check that the arguments are lists or Numpy arrays
    array = np.asarray(array, dtype=np.complex)
    kernel = np.asarray(kernel, dtype=np.complex)

    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise Exception('array and kernel have differing number of'
                        'dimensions')

    # store the dtype for conversion back later
    array_dtype = array.dtype
    # turn the arrays into 'complex' arrays
    if array.dtype.kind != 'c':
        array = array.astype(np.complex)
    if kernel.dtype.kind != 'c':
        kernel = kernel.astype(np.complex)

    # mask catching - masks must be turned into NaNs for use later
    if np.ma.is_masked(array):
        mask = array.mask
        array = np.array(array)
        array[mask] = np.nan
    if np.ma.is_masked(kernel):
        mask = kernel.mask
        kernel = np.array(kernel)
        kernel[mask] = np.nan

    # replace fftn if has_fftw so that nthreads can be passed
    global fftn, ifftn
    if has_fftw and not use_numpy_fft:
        def fftn(*args, **kwargs):
            return fftwn(*args, nthreads=nthreads, **kwargs)

        def ifftn(*args, **kwargs):
            return ifftwn(*args, nthreads=nthreads, **kwargs)
    elif use_numpy_fft:
        fftn = np.fft.fftn
        ifftn = np.fft.ifftn


    # NAN catching
    nanmaskarray = (array != array)
    array[nanmaskarray] = 0
    nanmaskkernel = (kernel != kernel)
    kernel[nanmaskkernel] = 0
    if ((nanmaskarray.sum() > 0 or nanmaskkernel.sum() > 0) and not interpolate_nan
            and not quiet):
        warnings.warn("NOT ignoring nan values even though they are present" +
                " (they are treated as 0)")

    if normalize_kernel is True:
        kernel = kernel / kernel.sum()
        kernel_is_normalized = True
    elif normalize_kernel:
        # try this.  If a function is not passed, the code will just crash... I
        # think type checking would be better but PEPs say otherwise...
        kernel = kernel / normalize_kernel(kernel)
        kernel_is_normalized = True
    else:
        if np.abs(kernel.sum() - 1) < 1e-8:
            kernel_is_normalized = True
        else:
            kernel_is_normalized = False


    if boundary is None:
        WARNING = ("The convolvend version of boundary=None is equivalent" +
                " to the convolve boundary='fill'.  There is no FFT " +
                " equivalent to convolve's zero-if-kernel-leaves-boundary" )
        warnings.warn(WARNING)
        psf_pad = True
    elif boundary == 'fill':
        # create a boundary region at least as large as the kernel
        psf_pad = True
    elif boundary == 'wrap':
        psf_pad = False
        fft_pad = False
        fill_value = 0 # force zero; it should not be used
    elif boundary == 'extend':
        raise NotImplementedError("The 'extend' option is not implemented " +
                "for fft-based convolution")

    arrayshape = array.shape
    kernshape = kernel.shape
    ndim = len(array.shape)
    if ndim != len(kernshape):
        raise ValueError("Image and kernel must " +
            "have same number of dimensions")
    # find ideal size (power of 2) for fft.
    # Can add shapes because they are tuples
    if fft_pad:
        if psf_pad:
            # add the dimensions and then take the max (bigger)
            fsize = 2**np.ceil(np.log2(
                np.max(np.array(arrayshape) + np.array(kernshape))))
        else:
            # add the shape lists (max of a list of length 4) (smaller)
            # also makes the shapes square
            fsize = 2**np.ceil(np.log2(np.max(arrayshape+kernshape)))
        newshape = np.array([fsize for ii in range(ndim)])
    else:
        if psf_pad:
            # just add the biggest dimensions
            newshape = np.array(arrayshape)+np.array(kernshape)
        else:
            newshape = np.array([np.max([imsh, kernsh])
                for imsh, kernsh in zip(arrayshape, kernshape)])


    # separate each dimension by the padding size...  this is to determine the
    # appropriate slice size to get back to the input dimensions
    arrayslices = []
    kernslices = []
    for ii, (newdimsize, arraydimsize, kerndimsize) in enumerate(zip(newshape, arrayshape, kernshape)):
        center = newdimsize - (newdimsize+1)//2
        arrayslices += [slice(center - arraydimsize//2,
            center + (arraydimsize+1)//2)]
        kernslices += [slice(center - kerndimsize//2,
            center + (kerndimsize+1)//2)]

    bigarray = np.ones(newshape, dtype=np.complex128) * fill_value
    bigkernel = np.zeros(newshape, dtype=np.complex128)
    bigarray[arrayslices] = array
    bigkernel[kernslices] = kernel
    arrayfft = fftn(bigarray)
    # need to shift the kernel so that, e.g., [0,0,1,0] -> [1,0,0,0] = unity
    kernfft = fftn(np.fft.ifftshift(bigkernel))
    fftmult = arrayfft * kernfft
    if (interpolate_nan or ignore_edge_zeros) and kernel_is_normalized:
        if ignore_edge_zeros:
            bigimwt = np.zeros(newshape, dtype=np.complex128)
        else:
            bigimwt = np.ones(newshape, dtype=np.complex128)
        bigimwt[arrayslices] = 1.0 - nanmaskarray * interpolate_nan
        wtfft = fftn(bigimwt)
        # I think this one HAS to be normalized (i.e., the weights can't be
        # computed with a non-normalized kernel)
        wtfftmult = wtfft * kernfft / kernel.sum()
        wtsm = ifftn(wtfftmult)
        # need to re-zero weights outside of the image (if it is padded, we
        # still don't weight those regions)
        bigimwt[arrayslices] = wtsm.real[arrayslices]
        # curiously, at the floating-point limit, can get slightly negative numbers
        # they break the min_wt=0 "flag" and must therefore be removed
        bigimwt[bigimwt < 0] = 0
    else:
        bigimwt = 1

    if np.isnan(fftmult).any():
        # this check should be unnecessary; call it an insanity check
        raise ValueError("Encountered NaNs in convolve.  This is disallowed.")

    # restore nans in original image (they were modified inplace earlier)
    # We don't have to worry about masked arrays - if input was masked, it was
    # copied
    array[nanmaskarray] = np.nan
    kernel[nanmaskkernel] = np.nan

    if return_fft:
        if fftshift:  # default on
            if crop:
                return np.fft.fftshift(fftmult)[arrayslices]
            else:
                return np.fft.fftshift(fftmult)
        else:
            return fftmult

    if interpolate_nan or ignore_edge_zeros:
        rifft = (ifftn(fftmult)) / bigimwt
        if not np.isscalar(bigimwt):
            rifft[bigimwt < min_wt] = np.nan
            if min_wt == 0.0:
                rifft[bigimwt == 0.0] = 0.0
    else:
        rifft = (ifftn(fftmult))

    if crop:
        result = rifft[arrayslices].real
        return result
    else:
        return rifft.real
