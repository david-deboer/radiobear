#! /usr/bin/env python
# script planet processing
import planet
import sys
import regrid
import shutil
import datetime


def splitFile(lineoutput='Scratch/specoutputline.dat'):
    fp = open(lineoutput, 'r')
    f = []
    T = []
    for line in fp:
        data = line.split()
        f.append(float(data[0]))
        del(data[0])
        ddd = []
        for d in data:
            ddd.append(float(d))
        T.append(ddd)
    fp.close()
    return f, T


oldFreqs = [8.46, 14.94, 22.46, 43.34]
newFreqs = [1.46, 2.44, 3.04, 3.44, 4.42, 5.45, 5.9, 6.45, 7.45, 8.58, 9.61, 10.38, 11.41,
            13.18, 14.21, 15.18, 16.21, 17.38, 22.45, 23.45, 24.45, 25.45]
newFreqsLo = [1.46, 2.44, 3.04, 3.44, 4.42, 5.45, 5.9]
newFreqsMi = [6.45, 7.45, 8.58, 9.61, 10.38, 11.41, 13.18]
newFreqsHi = [14.21, 15.18, 16.21, 17.38, 22.45, 23.45, 24.45, 25.45]

freqs = oldFreqs
freqs = 10.
print 'Reading input file ', sys.argv[1]

j = planet.planet('jupiter')
gasFile = 1
cloudFile = 1
if len(sys.argv) > 2:
    cloudFile = 2
j.atm.readGas(sys.argv[gasFile])
j.atm.readCloud(sys.argv[cloudFile])

j.atm.tweakAtm()
j.atm.computeProp()
regrid.regrid(j.atm, regridType=j.atm.config.regridType, Pmin=j.atm.config.pmin, Pmax=j.atm.config.pmax)
j.atm.nAtm = len(j.atm.gas[0])

b = '0.0,0.2,0.4,0.6,0.7,0.8,0.9<90'
j.run(freqs, b=b, outputType='batch')
nfreq, nTB = splitFile()
shutil.copy('specoutputline.dat', 'limb.dat')

j.run('reuse', 'disc', outputType='batch')
dfreq, dTB = splitFile()

fp = open('runpResults.dat', 'a')
fn = sys.argv[1].split('/')
s = fn[-1] + '\t'
fp.write(s)
s = ''
for i, f in enumerate(nfreq):
    if f != dfreq[i]:
        print 'Error', f, dfreq[i]
    s += '%.3f(' % (f)
    for T in nTB[i]:
        s += '%.3f,' % (T)
    for T in dTB[i]:
        s += '%.3f,' % (T)
    s = s.strip(',')
    s += ') '
s = s.strip()
s += '\n'
fp.write(s)
fp.close()

print datetime.datetime.now()
