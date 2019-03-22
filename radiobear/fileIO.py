from __future__ import print_function, absolute_import, division
import os.path
import numpy as np
import copy
from . import utils
import six


class FileIO(object):
    def __init__(self, directory):
        self.directory = directory
        self.header = {}

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

    def writeSpectrum(self, fp, freqs, freqUnit, b, Tb, xaxis='frequency', lineoutput='Scratch/specoutputline.dat'):
        fp_lineoutput = open(lineoutput, 'w')
        if xaxis.lower() == 'frequency':
            s = '# {}  K@b  \t'.format(freqUnit)
        elif xaxis.lower() == 'wavelength':
            s = '# cm   K@b  \t'
        else:
            s = '# {}    cm    K@b  \t'.format(freqUnit)
        for i, bv in enumerate(b):
            s += '({:5.3f},{:5.3f})\t'.format(bv[0], bv[1])
        s = s.strip('\t') + '\n'
        fp.write(s)
        for i, f in enumerate(freqs):
            wlcm = 100.0 * utils.speedOfLight / (f * utils.Units[freqUnit])
            if xaxis == 'frequency':
                s = '{:.2f}\t  '.format(f)
            elif xaxis == 'wavelength':
                s = '{:.4f}\t  '.format(wlcm)
            else:
                s = '{:.2f}     {:.4f} \t '.format(f, wlcm)
            for j in range(len(b)):
                s += '  {:9.4f}  \t'.format(Tb[j][i])
            s = s.strip() + '\n'
            fp_lineoutput.write(s)
            fp.write(s)
        fp_lineoutput.close()

    def writeProfile(self, fp, freqs, freqUnit, b, Tb, xaxis='frequency'):
        if xaxis == 'frequency':
            s = '# b  K@{} \t'.format(freqUnit)
        elif xaxis == 'wavelength':
            s = '# b  K@cm  \t'
        else:
            s = '# b  K@{},cm  \t'.format(freqUnit)
        for i, fv in enumerate(freqs):
            wlcm = 100.0 * utils.speedOfLight / (fv * utils.Units[freqUnit])
            if xaxis == 'frequency':
                s += '  {:9.4f}   \t'.format(fv)
            elif xaxis == 'wavelength':
                s += '  {:.4f}   \t'.format(wlcm)
            else:
                s += ' {:.2f},{:.4f}\t'.format(fv, wlcm)
        s = s.strip() + '\n'
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
            sfile = six.moves.input('File numbers: ')
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

    def read(self, fn=None, file_type='spectrum', directory=None):
        """reads brightness temperature file(s):
           fn = filename to read (but then ignores directory) | '?', '.' or None | integer [None]
           directory = subdirectory for data (not used if filename given) ['Output']"""

        if directory is None:
            directory = self.directory
        if fn is None:
            try_files = self.flist(fn, directory)
        else:
            if isinstance(fn, six.string_types):
                try_files = fn.split(',')
                try_files = [os.path.join(directory, x) for x in try_files]

        self.files = []
        self.TB = {}
        self.freqs = {}
        self.b = {}
        self.x = {}
        self.y = {}
        # ## Read in two passes:  first header
        for i, filename in enumerate(try_files):
            try:
                fp = open(filename, 'r')
            except IOError:
                print(filename + " not found - removing from list")
                continue
            if file_type.lower() != 'all' and file_type.lower() not in filename.lower():
                continue
            print("\tReading " + filename)
            self.files.append(filename)

            # ## Get past any header and get first line
            headerText = []
            for line in fp:
                if line[0] == '#':
                    headerText.append(line)
            fp.close()
            self.header[filename] = self.parseHeader(headerText)

        validData = False
        for filename in self.files:
            freqs = []  # frequencies
            TB = []   # data
            b = []      # b-vectors (resolution for images)
            # ## Now we have valid files and the header, now read in ftype
            with open(filename, 'r') as fp:
                # Get past header
                for line in fp:
                    if line[0] == '#':
                        if 'K@' in line:
                            label_line = copy.copy(line)
                        continue
                    else:
                        break
                types_in_file = [x.lower() for x in self.header[filename]['outType']]
                is_spectrum = 'spectrum' in types_in_file
                is_profile = 'profile' in types_in_file
                is_image = 'image' in types_in_file
                if is_image:
                    imRow = imCol = 0
                    b = self.resolution
                    bmag = b
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
                        TB.append(vdat)
                    self.header[filename]['img_filename'] = filename
                    if 'img_size' not in self.header[filename].keys():
                        self.header[filename]['img_size'] = [imRow, imCol]
                    else:
                        print(self.header[filename]['img_size'])
                        print('should equal ', imRow, imCol)
                    self.TB[filename] = np.array(TB)
                    xyextents = [-self.resolution * len(self.TB) / 2.0, self.resolution * len(self.TB) / 2.0,
                                 -self.resolution * len(self.TB) / 2.0, self.resolution * len(self.TB) / 2.0]
                    x = []
                    y = []
                    for i in range(len(self.TB)):
                        x.append(xyextents[0] + i * resolution)
                        y.append(xyextents[2] + i * resolution)
                    self.x[filename] = np.array(x)
                    self.y[filename] = np.array(y)
                    self.TB[filename] = np.array(TB)
                elif is_spectrum or is_profile:
                    validData = True
                    labels = label_line.split()
                    xlabel = labels[1]
                    ylabel = labels[2].split('@')[1]
                    del(labels[0:3])
                    curveLabels = labels
                    if is_spectrum:
                        btmp = []
                        print('b = ', end='')
                        for bb in labels:
                            print(' ' + bb, end='')
                            bbb = bb.split('(')[1].strip(')').split(',')
                            bbbb = [float(bbb[0]), float(bbb[1])]
                            btmp.append(bbbb)
                        b = np.array(btmp)
                        print('')
                    elif is_profile:
                        print('Freq = ', end='')
                        freq = []
                        for f in labels:
                            try:
                                ff = float(f)
                            except ValueError:
                                print(line)
                                continue
                            freqs.append(ff)
                            print('{:.3f} {}'.format(ff, 'GHz_hardcoded'), end='')
                        freqs = np.array(freqs)
                    for line in fp:
                        if line[0] == '#':
                            continue
                        elif len(line) > 2:
                            dat = line.split()
                            for i in range(len(dat)):
                                try:
                                    dat[i] = float(dat[i])
                                except ValueError:
                                    validData = False
                            if is_spectrum:
                                freqs.append(dat[0])
                            elif is_profile:
                                print("NEED TO READ B'S")
                            del(dat[0])
                            TB.append(dat)
                if validData:
                    self.TB[filename] = np.array(TB).transpose()
                    self.freqs[filename] = np.array(freqs)
                    self.b[filename] = np.array(b)

    def parseHeader(self, headerText):
        """Parses the pyPlanet image header"""
        _header = {}
        for hdr in headerText:
            hdr = hdr.strip('#').strip()
            updateKey = False
            if ':' in hdr:
                updateKey = True
                hdrkey = hdr.split(':')[0]
                hdr = ':'.join(hdr.split(':')[1:])
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
                if hdrkey not in _header.keys():
                    updateKey = True
                    h = [None]
            if updateKey:
                _header[hdrkey] = h
                self.__dict__[hdrkey] = h
        # ## set any header-derived values
        if 'res' in _header.keys():
            self.resolution = _header['res'][0]   # keep this for backward compatibility
        if 'freqs' in _header.keys():
            self.freq = _header['freqs'][0]
            _header['band'] = utils.getRFband(self.freq, 'GHz')
        return _header
