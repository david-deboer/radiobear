from __future__ import print_function, absolute_import, division
import os.path
import numpy as np
from . import utils
import six


class FileIO(object):
    def __init__(self, outputType, output_directory):
        self.outputType = outputType
        self.directory = output_directory

    def write(self, outputFile, outType, freqs, freqUnit, b, Tb, header):
        with open(outputFile, 'w') as fp:
            self.writeHeader(header, fp)
            if outType.lower() == 'image':
                self.writeImage(fp, Tb)
            elif outType.lower() == 'spectrum':
                self.writeSpectrum(fp, freqs, freqUnit, b, Tb)
            elif outType.lower() == 'profile':
                self.writeProfile(fp, freqs, freqUnit, b, Tb)
            else:
                print("Invalid output type: {}".format(outType))
        return outputFile

    def writeSpectrum(self, fp, freqs, freqUnit, b, Tb, lineoutput='Scratch/specoutputline.dat'):
        fp_lineoutput = open(lineoutput, 'w')
        if self.outputType.lower() == 'frequency':
            s = '# {}  K@b  \t'.format(freqUnit)
        elif self.outputType.lower() == 'wavelength':
            s = '# cm   K@b  \t'
        else:
            s = '# {}    cm    K@b  \t'.format(freqUnit)
        for i, bv in enumerate(b):
            s += '({:5.3f},{:5.3f})\t'.format(bv[0], bv[1])
        s = s.strip('\t') + '\n'
        fp.write(s)
        for i, f in enumerate(freqs):
            wlcm = 100.0 * utils.speedOfLight / (f * utils.Units[freqUnit])
            if self.outputType == 'frequency':
                s = '{:.2f}\t  '.format(f)
            elif self.outputType == 'wavelength':
                s = '{:.4f}\t  '.format(wlcm)
            else:
                s = '{:.2f}     {:.4f} \t '.format(f, wlcm)
            for j in range(len(b)):
                s += '  {:7.2f}  \t'.format(Tb[j][i])
            s = s.strip() + '\n'
            fp_lineoutput.write(s)
            fp.write(s)
        fp_lineoutput.close()

    def writeProfile(self, fp, freqs, freqUnit, b, Tb):
        if self.outputType == 'frequency':
            s = '# b  K@{} \t'.format(freqUnit)
        elif self.outputType == 'wavelength':
            s = '# b  K@cm  \t'
        else:
            s = '# b  K@{},cm  \t'.format(freqUnit)
        for i, fv in enumerate(freqs):
            wlcm = 100.0 * utils.speedOfLight / (fv * utils.Units[freqUnit])
            if self.outputType == 'frequency':
                s += '  {:9.4f}   \t'.format(fv)
            elif self.outputType == 'wavelength':
                s += '  {:.4f}   \t'.format(wlcm)
            else:
                s += ' {:.2f},{:.4f}\t'.format(fv, wlcm)
        s.strip('\t') + '\n'
        fp.write(s)
        bs = []
        for i, bv in enumerate(b):
            s = '{:5.3f} {:5.3f}\t'.format(bv[0], bv[1])
            bs.append(np.sqrt(bv[0]**2 + bv[1]**2))
            for j in range(len(freqs)):
                s += ' {:7.2f}\t '.format(Tb[i][j])
            s = s.strip() + '\n'
            fp.write(s)

    def writeImage(self, fp, Tb):
        for data in Tb:
            s = ''
            for d in data:
                s += '{:7.2f}\t'.format(d)
            s = s.strip() + '\n'
            fp.write(s)

    def writeHeader(self, header, fp):
        alpha_header = sorted(header.keys())
        for hdr in alpha_header:
            fp.write(header[hdr])

    def flist(self, fd=None, directory=None, tag='dat'):
        """This generates the list of filenames to be opened - doesn't check for existence"""
        if directory is not None:
            usedir = directory
        else:
            usedir = self.directory
        file_list = utils.ls(directory=usedir, tag=tag, show=False, returnList=True)

        ifile = []
        files = []
        if isinstance(fd, six.integer_types):
            ifile = [fd]
        elif fd is None:
            for i, fn in enumerate(file_list):
                print('{}  -  {}'.format(i, fn))
            sfile = raw_input('File numbers: ')
            split_on = None
            if '-' in sfile:
                sfile = sfile.split('-')
                ifile = range(int(sfile[0]), int(sfile[1]) + 1)
            else:
                if ',' in sfile:
                    split_on = ','
                ifile = [int(x) for x in sfile.split(split_on)]
        elif type(fd) == list and isinstance(fd[0], six.integer_types):
            ifile = fd
        elif type(fd) == list and isinstance(fd[0], six.string_types):
            files = fd
        else:
            files = [fd]

        if bool(len(ifile)):
            for i, fn in enumerate(file_list):
                if i in ifile:
                    files.append(fn)

        return files

    def read(self, fd=None, directory=None):
        """reads brightness temperature file(s):
           fn = filename to read (but then ignores directory) | '?', '.' or None | integer [None]
           directory = subdirectory for data (not used if filename given) ['Output']"""

        ftypes = ['image', 'profile', 'spectrum']
        if fd is None:
            if directory is not None:
                usedir = directory
            else:
                usedir = self.directory
            try_files = self.flist(fd, usedir)
        else:
            try_files = fd

        # ## Read in two passes:  first header
        self.ftype = None
        headerText = []
        self.files = []
        self.ftype = 'Unknown'
        for i, filename in enumerate(try_files):
            try:
                fp = open(filename, 'r')
            except IOError:
                print(filename + " not found - removing from list")
                continue
            for ff in ftypes:
                if ff in filename:
                    self.ftype = ff
                    break
            if not i:
                overallType = self.ftype
                print('File type:  ' + overallType)
            else:
                if self.ftype != overallType:
                    print(self.ftype + ' not of type ' + overallType + ' - removing from list')
                    continue

            print("\tReading " + filename)
            self.files.append(filename)

            # ## Get past any header and get first line
            line = '# file:  ' + filename + '\n'
            for line in fp:
                if line[0] == '#':
                    headerText.append(line)
            fp.close()
        self.parseHeader(headerText)

        freqs = []  # frequencies
        wavel = []  # wavelengths
        data = []   # data
        b = []      # b-vectors (resolution for images)
        bmag = []   # b-mag (resolution for images)
        # ## Now we have valid files and the header, now read in ftype
        if self.ftype == 'image':
            imRow = imCol = 0
            b = self.resolution
            bmag = b
            for filename in self.files:
                fp = open(filename, 'r')
                for line in fp:
                    if line[0] == '#':
                        continue
                    imRow += 1
                    vdat = []
                    sdat = line.split()
                    if len(sdat) < 2:
                        continue
                    imCol = 0
                    for v in sdat:
                        imCol += 1
                        vdat.append(float(v))
                    data.append(vdat)
                fp.close()
            self.header['img_filename'] = filename
            if 'img_size' not in self.header.keys():
                self.header['img_size'] = [imRow, imCol]
            else:
                print(self.header['img_size'])
                print('should equal ', imRow, imCol)
            self.data = np.array(data)
            self.xyextents = [-self.resolution * len(self.data) / 2.0, self.resolution * len(self.data) / 2.0,
                              -self.resolution * len(self.data) / 2.0, self.resolution * len(self.data) / 2.0]
            self.x = []
            self.y = []
            for i in range(len(self.data)):
                self.x.append(self.xyextents[0] + i * self.resolution)
                self.y.append(self.xyextents[2] + i * self.resolution)
            self.x = np.array(self.x)
            self.y = np.array(self.y)
        elif self.ftype == 'spectrum' or self.ftype == 'profile':
            # indata = []
            n = 0
            validData = True
            for filename in self.files:
                fp = open(filename, 'r')
                for line in fp:
                    if line[0] == '#' and 'K@' in line:
                        labels = line.split()
                        xlabel = labels[1]
                        ylabel = labels[2].split('@')[1]
                        del(labels[0:3])
                        curveLabels = labels
                        if self.ftype == 'spectrum':
                            print('b = ', end='')
                            for bb in labels:
                                print(' ' + bb, end='')
                                bbb = bb.split('(')[1].strip(')').split(',')
                                bbbb = [float(bbb[0]), float(bbb[1])]
                                b.append(bbbb)
                                bmag.append(np.sqrt(bbbb[0]**2 + bbbb[1]**2))
                            b = np.array(b)
                            bmag = np.array(bmag)
                            print('')
                        elif self.ftype == 'profile':
                            print('Freq = ', end='')
                            for f in labels:
                                try:
                                    ff = float(f)
                                except ValueError:
                                    print(line)
                                    continue
                                freqs.append(ff)
                                print('{:.3f} {}'.format(ff, 'GHz_hardcoded'), end='')
                            freqs = np.array(freqs)
                    elif line[0] == '#':
                        continue
                    elif len(line) > 2:
                        dat = line.split()
                        for i in range(len(dat)):
                            try:
                                dat[i] = float(dat[i])
                            except ValueError:
                                dat[i] = dat[i]
                                validData = False
                        if self.ftype == 'spectrum':
                            freqs.append(dat[0])
                            wavel.append((utils.speedOfLight / 1.0E7) / dat[0])
                        elif self.ftype == 'profile':
                            print("NEED TO READ B'S")
                        del(dat[0])
                        data.append(dat)
                fp.close()
        if validData:
            self.data = np.array(data)
        self.freqs = np.array(freqs)
        self.wavel = np.array(wavel)
        self.b = np.array(b)
        self.bmag = np.array(bmag)

    def parseHeader(self, headerText):
        """Parses the pyPlanet image header"""
        for hdr in headerText:
            hdr = hdr.strip('#').strip()
            print(hdr)
            updateKey = False
            if ':' in hdr:
                updateKey = True
                hdrkey = hdr.split(':')[0].split()[0]
                hdr = hdr.split(':')[1]
                h = []
                hdr = hdr.split()
                for dat in hdr:
                    dat = dat.strip('[').strip(']').strip(',')
                    try:
                        h.append(float(dat))
                    except ValueError:
                        h.append(dat)
            else:
                hdrkey = hdr.split()[0]
                if hdrkey not in self.header.keys():
                    updateKey = True
                    h = [None]
            if updateKey:
                self.header[hdrkey] = h
                self.__dict__[hdrkey] = h
        # ## set any header-derived values
        if 'res' in self.header.keys():
            self.resolution = self.header['res'][0]   # keep this for backward compatibility
        if 'freqs' in self.header.keys():
            self.freq = self.header['freqs'][0]
            self.header['band'] = utils.getRFband(self.freq, 'GHz')
