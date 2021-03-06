import ph3_jh
import matplotlib.pyplot as plt
import numpy as np

f = np.arange(1.0, 60.0, 0.001)

P = 0.1
T = 300.0
X_ph3 = 3.3E-4
X_he = 0.1
X_h2 = 1.0 - (X_ph3+X_he)
X = [X_h2, X_he, X_ph3]
P_dict = {'H2': 0, 'HE': 1, 'PH3': 2}
a = ph3_jh.alpha(f, T, P, X, P_dict, None, truncate_strength=None, truncate_freq=None)
plt.semilogy(f, a)
