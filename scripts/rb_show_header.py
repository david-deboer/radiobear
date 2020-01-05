#! /usr/bin/env python
import argparse
from radiobear import fileIO

ap = argparse.ArgumentParser()

ap.add_argument('file', help="File from which to show header.")
ap.add_argument('-d', '--directory', help="Directory in which to look.", default="Output")
ap.add_argument('--file-type', dest='file_type', help="file-type to show.", default='all')
args = ap.parse_args()

fio = fileIO.FileIO(directory=args.directory)
fio.read(fn=args.file, file_type=args.file_type)

for f in fio.files:
    for key in fio.header[f].keys():
        print("{}:  {}".format(key, fio.header[f][key]))
