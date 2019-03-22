#! /usr/bin/env python
import os
import argparse

ap = argparse.ArgumentParser()

ap.add_argument('tag', help="Tag to delete in directory")
ap.add_argument('-d', '--directory', help="Directory in which to delete tagged files", default=None)
args = ap.parse_args()

if args.directory is None:
    if args.tag == 'dat':
        args.directory = 'Output'
    elif args.tag == 'log':
        args.directory = 'Logs'

print("This will delete all {} files in {}".format(args.tag, args.directory))
are_you_sure = six.moves.input("Are you sure (y/n)?")
if are_you_sure == 'y':
    files_list = os.listdir(args.directory)
    for fa in files_list:
        tag = fa.split('.')[-1]
        if tag == args.tag:
            os.remove(os.path.join(args.directory, fa))
