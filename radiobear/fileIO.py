# -*- mode: python; coding: utf-8 -*-
"""FileIO."""
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.

import numpy as np
from . import utils
from . import data_handling


class FileIO(object):
    """File input/output class."""

    def __init__(self, directory='Output', scratch_spec='Scratch/specoutputline.dat'):
        """Initialize."""
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
        """Write spectrum file."""
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
        """Write profile file."""
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
        """Write image file."""
        for data in self.data.Tb:
            s = ''
            for d in data:
                s += '{:7.2f}\t'.format(d)
            s = s.strip() + '\n'
            fp.write(s)

    def writeHeader(self, header, fp):
        """Write header."""
        alpha_header = sorted(header.keys())
        for hdr in alpha_header:
            fp.write(header[hdr] + '\n')

    def flist(self, files=None, tag=None):
        """
        Generate the list of filenames to be opened.

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
        if isinstance(files, int):
            ifile = [files]
        elif isinstance(files, str):
            file_list = files.split(',')
            ifile = range(len(file_list))
        elif isinstance(files, list):
            if isinstance(files[0], int):
                ifile = files[:]
            elif isinstance(files[0], str):
                file_list = files[:]
                ifile = range(len(files))
        elif files is None:
            for i, fn in enumerate(file_list):
                print('{}  -  {}'.format(i, fn))
            sfile = input('File numbers (n; n1-n2; n1,n2,...; all): ')
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

    def show(self, property='all'):
        """Show properties."""
        data_dict = vars(self.data)
        if isinstance(property, str):
            if property == 'all':
                keys = list(data_dict.keys())
            else:
                keys = [property]
        else:
            keys = property
        for k in keys:
            if k == 'header':
                self.data.show_header()
            elif k == 'log':
                self.data.show_log()
            else:
                print("<<<{}>>>:  {}".format(k, data_dict[k]))

    def read(self, fn=None, tag='dat', file_type='spectrum'):
        """
        Read brightness temperature file(s).

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
        print('\n')
        for i, filename in enumerate(try_files):
            try:
                fp = open(filename, 'r')
            except IOError:
                print(filename + " not found - removing from list")
                continue
            if file_type.lower() != 'all' and file_type.lower() not in filename.lower():
                continue
            print("Reading " + filename)
            self.files.append(filename)

            # ## Get past any header and get first line
            header_text = []
            for line in fp:
                if line[0] == '#':
                    header_text.append(line)
            fp.close()
            self.data[filename] = data_handling.Data()
            self.data[filename].header = self.read_header(header_text)
            self.set_header_attr(filename)
            self.data[filename].show(include=['header', 'start', 'stop', 'logfile'])

        for filename in self.files:
            self.data[filename].f = []  # frequencies
            self.data[filename].Tb = []  # data
            self.data[filename].b = []  # b-vectors (resolution for images)
            type_in_file = self.data[filename].header['type'].lower()
            is_type = {'spectrum': 'spectrum' in type_in_file,
                       'profile': 'profile' in type_in_file,
                       'image': 'image' in type_in_file}
            if is_type['image']:
                valid_data = self._process_image(filename)
            else:
                valid_data = self._process_other(filename, is_type)
            if valid_data:
                self.data[filename].Tb = np.array(self.data[filename].Tb).transpose()
                self.data[filename].f = np.array(self.data[filename].f)
                self.data[filename].b = np.array(self.data[filename].b)

    def _process_other(self, filename, is_type):
        label_line_processed = False
        with open(filename, 'r') as fp:
            xxx = []
            for line in fp:
                print(line)
                if line[0] == '#':
                    if 'K@' in line:
                        label_line = line[1:].split()
                    continue
                if not label_line_processed:
                    label_line_processed = True
                    labels = label_line[2:]
                    # Process label_line
                    if is_type['spectrum']:
                        nl = []
                        for plbl in labels:
                            ltmp = plbl.replace(')', ' ').replace('(', '').strip()
                            nl.append([float(x) for x in ltmp.split(',')])
                        self.data[filename].b = nl
                        print('b = ', self.data[filename].b)
                    elif is_type['profile']:
                        self.data[filename].f = [float(x) for x in labels]
                        print('Freq = ', self.data[filename].f)

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

    def _process_image(self, fp, filename, resolution):
        print("NOT IMPLEMENTED YET")
        return False
        imRow = imCol = 0
        self.data[filename].b = self.resolution
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
        xyextents = [-self.resolution * len(self.data[filename].Tb) / 2.0, self.resolution *
                     len(self.data[filename].Tb) / 2.0,
                     -self.resolution * len(self.data[filename].Tb) / 2.0, self.resolution *
                     len(self.data[filename].Tb) / 2.0]
        x = []
        y = []
        for i in range(len(self.TB)):
            x.append(xyextents[0] + i * resolution)
            y.append(xyextents[2] + i * resolution)
        self.data[filename].x = np.array(x)
        self.data[filename].y = np.array(y)

    def set_header_attr(self, filename):
        """Set header attribute."""
        for key in self.starred_header_keys:
            if key in self.data[filename].allowed_parameters:
                self.data[filename].set(key, self.data[filename].header[key])

    def read_header(self, header_text):
        """Read the radiobear header."""
        _header = {}
        self.starred_header_keys = []
        for hdr in header_text:
            hdrproc = hdr.strip('#').strip()
            if 'K@' in hdr:
                _header['label-line'] = hdrproc
                continue
            split_on = ':' if ':' in hdr else None
            hdrkey = hdrproc.split(split_on)[0]
            if split_on is not None:
                hdrproc = hdrproc[(hdrproc.find(split_on) + 1):].strip()
            if hdrkey in _header.keys():
                print("'{}' already included. Pushed to 'other_'".format(hdr, hdrkey))
                hdrkey = 'other_{}'.format(hdrkey)
            if hdrkey.strip().startswith('*'):
                hdrkey = hdrkey.strip().strip('*').strip()
                self.starred_header_keys.append(hdrkey)
            _header[hdrkey] = hdrproc
        # ## set any header-derived values
        if 'res' in _header.keys():
            self.resolution = _header['res'][0]   # keep this for backward compatibility
        return _header
