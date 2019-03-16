import clouds_idp as cl
import matplotlib.pyplot as plt
import cmath
import math
import numpy as np

freqs = []
for i in range(10000):
    freqs.append(float(i + 1) / 30.0)


def e(Ts=[100.0, 150.0, 200.0, 250.0, 300.0, 350.0]):
    FR = [1.0E8, 3.0E8, 1.0E9, 2.0E9, 3.0E9, 5.0E9, 1.0E10, 3.0E10, 1.0E11]
    ff = np.array(FR) / 1.0E9
    EIMAG = [8.0E-3, 1.5E-3, 8.0E-4, 1.0E-3, 1.2E-3, 1.5E-3, 3.0E-3, 8.0E-3, 2.0E-2]
    ei = np.array(EIMAG)
    plt.semilogx(ff, ei, 'o')
    for T in Ts:
        d1 = []
        d2 = []
        for f in freqs:
            n = cl.refractiveIndex(f, T)
            e = n**2
            d1.append(e.real)
            d2.append(-e.imag)
        plt.semilogx(freqs, d1, label='d1 @ %.0f K' % (T))
        plt.semilogx(freqs, d2, label='d2 @ %.0f K' % (T))


P = 1.0
X = [1E-6, 1E-6]
cdict = {'H2O': 0, 'SOLN': 1}

other = {'water': 1}
T = 273.1
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))
T = 283.0
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))
T = 293.0
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))

other = {'ice': 1}
T = 272.9
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))
T = 268.0
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))
T = 263.0
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))
T = 253.0
a = cl.alpha(freqs, T, P, X, cdict, other)
plt.loglog(freqs, a, label='%f' % (T))
