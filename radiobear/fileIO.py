from __future__ import print_function, absolute_import, division
import os.path
import numpy as np
import copy
from . import utils
from . import data_handling
import six


class FileIO(object):
    def __init__(self, directory='Output', scratch_spec='Scratch/specoutputline.dat'):
        self.directory = directory
        self.scratch_spec = scratch_spec

    def write(self, output_file, data, xaxis='frequency'):
        """
        Write data to 'output_file'.

        Parameters:
        ------------
        output_file:  name of output file
        data:  class data_handling.Data containing everything to write.
        xaxis:  'frequency' or 'wavelength'
        """
        self.data = data
        xaxis = xaxis.lower()
        with open(output_file, 'w') as fp:
            self.writeHeader(data.header, fp)
            if data.type.lower() == 'image':
                self.writeImage(fp)
            elif data.type.lower() == 'spectrum':
                self.writeSpectrum(fp, xaxis=xaxis)
            elif data.type.lower() == 'profile':
                self.writeProfile(fp, xaxis=xaxis)
            else:
                print("Invalid output type: {}".format(data.type))
        return output_file

    def writeSpectrum(self, fp, xaxis):
        fp_lineoutput = open(self.scratch_spec, 'w')
        if xaxis.startswith('wave'):
            s = '# cm   K@b  \t'
        else:
            s = '# {}  K@b  \t'.format(self.data.freqUnit)
        if utils.b_type(self.data.b) == 'disc':
            s += 'disc'
        else:
            for i, bv in enumerate(self.data.b):
                s += '({:5.3f},{:5.3f})\t'.format(bv[0], bv[1])
        s = s.strip('\t') + '\n'
        fp.write(s)
        for i, f in enumerate(self.data.f):
            wlcm = 100.0 * utils.speedOfLight / (f * utils.Units[self.data.freqUnit])
            if xaxis.startswith('wave'):
                s = '{:.4f}\t  '.format(wlcm)
            else:
                s = '{:.2f}\t  '.format(f)
            for j in range(len(self.data.b)):
                s += '  {:9.4f}  \t'.format(self.data.Tb[j][i])
            s = s.strip() + '\n'
            fp_lineoutput.write(s)
            fp.write(s)
        fp_lineoutput.close()

    def writeProfile(self, fp, xaxis):
        if xaxis.startswith('wave'):
            s = '# b  K@cm  \t'
        else:
            s = '# b  K@{} \t'.format(self.data.freqUnit)
        for i, f in enumerate(self.data.f):
            wlcm = 100.0 * utils.speedOfLight / (f * utils.Units[self.data.freqUnit])
            if xaxis.startswith('wave'):
                s += '  {:.4f}   \t'.format(wlcm)
            else:
                s += '  {:9.4f}   \t'.format(f)
        s = s.strip() + '\n'
        fp.write(s)
        for i, bv in enumerate(self.data.b):
            if utils.b_type(self.data.b) == 'disc':
                s = 'disc'
            else:
                s = '{:5.3f} {:5.3f}\t'.format(bv[0], bv[1])
            for j in range(len(self.data.f)):
                s += ' {:7.2f}\t '.format(self.data.Tb[i][j])
            s = s.strip() + '\n'
            fp.write(s)

    def writeImage(self, fp):
        for data in self.data.Tb:
            s = ''
            for d in data:
                s += '{:7.2f}\t'.format(d)
            s = s.strip() + '\n'
            fp.write(s)

    def writeHeader(self, header, fp):
        alpha_header = sorted(header.keys())
        for hdr in alpha_header:
            fp.write(header[hdr])

    def flist(self, files=None, tag=None):
        """
        This generates the list of filenames to be opened.

        Parameters:
        ------------
        files:  file(s):  None, <str>, <int>, <list>
        tag:  tag to filter on:  <str> or None
        """
        usedir = self.directory
        file_list = utils.ls(directory=usedir, tag=tag, show=False, returnList=True)
        if not len(file_list):
            print("No files found in {}".format(usedir))
            return []
        file_list = sorted(file_list)

        ifile = []  # list of integers of files to use within file_list
        if isinstance(files, six.integer_types):
            ifile = [files]
        elif isinstance(files, six.string_types):
            file_list = files.split(',')
            ifile = range(len(file_list))
        elif isinstance(files, list):
            if isinstance(files[0], six.integer_types):
                ifile = files[:]
            elif isinstance(files[0], six.string_types):
                file_list = files[:]
                ifile = range(len(files))
        elif files is None:
            for i, fn in enumerate(file_list):
                print('{}  -  {}'.format(i, fn))
            sfile = six.moves.input('File numbers (n; n1-n2; n1,n2,...; all): ')
            if sfile.lower() == 'all':
                ifile = range(len(file_list))
            elif '-' in sfile:
                sfile = sfile.split('-')
                ifile = range(int(sfile[0]), int(sfile[1]) + 1)
            else:
                ifile = [int(x) for x in sfile.split(',')]

        ret_files = []
        for i in ifile:
            ret_files.append(file_list[i])
        return ret_files

    def read(self, fn=None, tag='dat', file_type='spectrum'):
        """
        Reads brightness temperature file(s):

        Parameters:
        ------------
        files:  file(s):  None, <str>, <int>, <list>
        tag:  tag to filter on:  <str> or None
        file_type:  'spectrum' or 'profile' etc
        """
        try_files = self.flist(fn, tag)
        self.files = []
        self.logs = []
        self.data = {}
        # ## Read in two passes:  first header and get files
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
            header_text = []
            for line in fp:
                if line[0] == '#':
                    header_text.append(line)
            fp.close()
            self.data[filename] = data_handling.Data()
            self.data[filename].header = self.parseHeader(header_text)
        for k, v in six.iteritems(self.data[filename].header):
            if k in self.data[filename].allowed_parameters:
                setattr(self.data[filename], k, v)

        for filename in self.files:
            self.data[filename].f = []  # frequencies
            self.data[filename].Tb = []  # data
            self.data[filename].b = []  # b-vectors (resolution for images)
            with open(filename, 'r') as fp:
                # Get past header and get label_line
                for line in fp:
                    if line[0] == '#':
                        if 'K@' in line:
                            label_line = copy.copy(line)
                        continue
                    else:
                        break
                types_in_file = [x.lower() for x in self.data[filename].header['data_type']]
                is_type = {'spectrum': 'spectrum' in types_in_file,
                           'profile': 'profile' in types_in_file,
                           'image': 'image' in types_in_file}
                if is_type['image']:
                    valid_data = self._process_image(fp, filename)
                else:
                    valid_data = self._process_other(fp, filename, is_type, label_line)
                if valid_data:
                    self.data[filename].Tb = np.array(self.data[filename].Tb).transpose()
                    self.data[filename].f = np.array(self.data[filename].f)
                    self.data[filename].b = np.array(self.data[filename].b)

    def _process_other(self, fp, filename, is_type, label_line):
        labels = label_line.split()
        xlabel = labels[1]
        ylabel = labels[2].split('@')[1]
        del(labels[0:3])
        # Process label_line
        if is_type['spectrum']:
            self.data[filename].b = labels
            print('b = ', self.data[filename].b)
            if not isinstance(self.data[filename].b[0], six.string_types):
                self.data[filename].b = [float(x) for x in self.data[filename].b]
        elif is_type['profile']:
            self.data[filename].f = labels.split()
            print('Freq = ', self.data[filename].f)
            if not isinstance(self.data[filename].f[0], six.string_types):
                self.data[filename].f = [float(x) for x in self.data[filename].f]
        # Process data
        xxx = []
        for line in fp:
            if line[0] == '#' or len(line) < 2:
                continue
            try:
                dat = [float(x) for x in line.split()]
            except ValueError:
                return False
            xxx.append(dat[0])
            del(dat[0])
            self.data[filename].Tb.append(dat)
        if is_type['spectrum']:
            self.data[filename].f = xxx
        else:
            self.data[filename].b = xxx
        return True

    def _process_image(self, fp, filename):
        print("NOT IMPLEMENTED YET")
        return False
        imRow = imCol = 0
        self.data[filename].b = self.resolution
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
            self.data[filename].Tb.append(vdat)
        self.data[filename].header['img_filename'] = filename
        if 'img_size' not in self.data[filename].header.keys():
            self.data[filename].header['img_size'] = [imRow, imCol]
        else:
            print(self.data[filename].header['img_size'])
            print('should equal ', imRow, imCol)
        self.data[filename].Tb = np.array(self.data[filename].Tb)
        xyextents = [-self.resolution * len(self.data[filename].Tb) / 2.0, self.resolution * len(self.data[filename].Tb) / 2.0,
                     -self.resolution * len(self.data[filename].Tb) / 2.0, self.resolution * len(self.data[filename].Tb) / 2.0]
        x = []
        y = []
        for i in range(len(self.TB)):
            x.append(xyextents[0] + i * resolution)
            y.append(xyextents[2] + i * resolution)
        self.data[filename].x = np.array(x)
        self.data[filename].y = np.array(y)

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
