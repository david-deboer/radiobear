#! /usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import argparse
from radiobear import utils

ap = argparse.ArgumentParser()
ap.add_argument('constituent', help='Name of constituent')
ap.add_argument('-x', '--mixing', help='Mixing ratio of constituent', default=1e-3)
ap.add_argument('-f', '--freqs', help='fstart,fstop,fstep in GHz', default='1,100,1')
ap.add_argument('--he', help='Mixing ratio of helium', default=0.1)
ap.add_argument('-p', '--pressure', help='Pressure in bars', default=1.0)
ap.add_argument('-t', '--temperature', help='Temperature in K', default=300.0)
args = ap.parse_args()
args.constituent = args.constituent.upper()
freqs = [float(x) for x in args.freqs.split(',')]

rbc_path = utils.rb_path('constituents/{}'.format(args.constituent.lower()))
if args.constituent == 'PH3':
    from radiobear.constituents.ph3 import ph3_jh
    alpha = ph3_jh.alpha
if args.constituent == 'NH3':
    from radiobear.constituents.nh3 import nh3_dbs_sjs
    alpha = nh3_dbs_sjs.alpha
if args.constituent == 'H2S':
    from radiobear.constituents.h2s import h2s_ddb
    alpha = h2s_ddb.alpha
if args.constituent == 'H2O':
    from radiobear.constituents.h2o import h2o_bk
    alpha = h2o_bk.alpha.alpha
if args.constituent == 'H2':
    from radiobear.constituents.h2 import h2_jj_ddb
    alpha = h2_jj_ddb.alpha
if args.constituent == 'CO':
    from radiobear.constituents.co import co_ddb
    alpha = co_ddb.alpha


f = np.arange(freqs[0], freqs[1], freqs[2])
P = float(args.pressure)
T = float(args.temperature)
X = float(args.mixing)
X_he = float(args.he)
X_h2 = 1.0 - (X - X_he)
X = [X_h2, X_he, X]
P_dict = {'H2': 0, 'HE': 1, args.constituent: 2}
a = alpha(f, T, P, X, P_dict, None, truncate_strength=None, truncate_freq=None, path=rbc_path)
plt.semilogy(f, a)
plt.show()
