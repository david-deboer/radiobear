import h2s_ddb as absorber
import matplotlib.pyplot as plt
import math
import numpy as np

f = []
fmin = 1.0
fmax = 500.0
fstep = 1.0
otherPar = []

f = np.arange(fmin, fmax, fstep)
P = 1.009
T = 216.4
X_absorber = 0.0095
X_he = 0.1347
X_h2 = 1.0 - (X_absorber + X_he)
X_partial = [X_h2, X_he, X_absorber]
P_dict = {'H2': 0, 'HE': 1, 'H2S': 2}
a = absorber.alpha(f, T, P, X_partial, P_dict, otherPar)
plt.semilogy(f, a)
plt.axis([1, 25, .01, 000])


plt.figure('vsP')
f = [0.5]
T = 200.0
Pvals = np.arange(1.0, 100.0, 1.0)
alpha = []
for P in Pvals:
    alpha.append(absorber.alpha(f, T, P, X_partial, P_dict, otherPar))
plt.plot(Pvals, alpha)

plt.figure('vsT')
f = [0.5]
P = 200.0
Tvals = np.arange(100.0, 1000.0, 10.0)
alpha = []
for T in Tvals:
    alpha.append(absorber.alpha(f, T, P, X_partial, P_dict, otherPar))
plt.plot(Tvals, alpha)
