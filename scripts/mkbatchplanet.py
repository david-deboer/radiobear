#! /usr/bin/env python
import os
import sys

print """Process is:
    make filenames.txt with input files
    'mkbatchplanet.py' with filenames in filenames.txt -- be sure to set range as arguments
        ==> writes 'batchplanet'
    set any parameters (use.py, freqs in runp.py, ...)
    'batchplanet'     does runp.py with filenames from filenames.txt (use screen -S sessionname batchplanet)
        ==> writes 'runpResults.dat'"""
print
print

path = '/Users/ddeboer/Documents/Projects/Planets/pyplanet/Jupiter/mike/set1'
path = '/indirect/o/ddeboer/pyplanet/Jupiter/mike/set1'
print path

names_are_set = False
try:
    filenames = open('filenames.txt','r')
except IOError:
    print 'filenames.txt not present'
else:
    names_are_set = True
    names = []
    for line in filenames:
        data = line.split()
        try:
            names.append(data[-1])
        except IndexError:
            print data+"  didn't work"

indices_are_set = False
if names_are_set:
    if len(sys.argv) == 1:
        print 'Default to use all files'
        number_to_use = 'all'
        indices_are_set = True
    elif len(sys.argv) == 2:
        try:
            start_at = 0
            number_to_use = int(sys.argv[1])
        except ValueError:
            number_to_use = 'all'
        indices_are_set = True
    elif len(sys.argv) == 3:
        try:
            start_at = int(sys.argv[1])
            number_to_use = int(sys.argv[2])
            indices_are_set = True
        except ValueError:
            indices_are_set = False
    else:
        indices_are_set = False
            
if indices_are_set:
    if number_to_use == 'all':
        start_at = 0
        number_to_use = len(names)
    end_at = start_at+number_to_use
    batchfile = open('batchplanet','w')
    use_names = names[start_at:end_at]
    s = "echo 'start' >> date.out\n"
    s+= "date >> date.out\n"
    batchfile.write(s)
    print "Writing 'batchplanet' with names for "+str(len(use_names))+" files (",
    print str(start_at)+" - "+str(end_at-1)+" )"
    for name in use_names:
        ffn = os.path.join(path,name)
        s = 'runp.py '+ffn+'\n'
        batchfile.write(s)
    s = "echo 'stop' >> date.out\n"
    s+= "date >> date.out\n"
    batchfile.write(s)
    os.chmod('batchplanet',0777)
else:
    print 'Did not run'
