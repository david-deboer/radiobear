#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('constituent', default='ph3')
args = ap.parse_args()
args.constituent = args.constituent.upper()

if args.constituent == 'PH3':
    from radiobear.constituents import ph3_jh
    alpha = ph3_jh.alpha
if args.constituent == 'NH3':
    from radiobear.constituents import nh3_dbs_sjs
    alpha = nh3_dbs_sjs.alpha
if args.constituent == 'H2S':
    from radiobear.constituents import h2s_ddb
    alpha = h2s_ddb.alpha
if args.constituent == 'H2O':
    from radiobear.constituents import h2o_bk
    alpha = h2o_bk.alpha.alpha
if args.constituent == 'H2':
    from radiobear.constituents import h2_jj_ddb
    alpha = h2_jj_ddb.alpha
if args.constituent == 'CO':
    from radiobear.constituents import co_ddb
    alpha = co_ddb.alpha


f = np.arange(1.0, 60.0, 0.001)

P = 0.1
T = 300.0
X_ph3 = 3.3E-4
X_he = 0.1
X_h2 = 1.0 - (X_ph3+X_he)
X = [X_h2, X_he, X_ph3]
P_dict = {'H2': 0, 'HE': 1, 'PH3': 2}
a = alpha(f, T, P, X, P_dict, None, truncate_strength=None, truncate_freq=None)
plt.semilogy(f, a)
