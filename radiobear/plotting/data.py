import matplotlib.pyplot as plt
import numpy as np
from radiobear import utils


# ##############################################################################################################
#                                          DATA RETURN PLOTTING (f, b, Tb)
# ##############################################################################################################

class plots:
    def __init__(self, data):
        self.data = data

    def show(self):
        plt.show()

    def profile(self):
        if utils.b_type(self.data.b) == 'disc':
            print("Disc average")
            return
        plt.figure('Profile')
        b = self.data.b.transpose()
        bvec = np.sqrt(b[0] * b[0] + b[1] * b[1])
        Tvec = np.transpose(self.data.Tb)
        for i, f in enumerate(self.data.f):
            plt.plot(bvec, Tvec[i], label="{:.2f} {}".format(f, self.data.freqUnit))
        plt.legend()
        plt.xlabel('b')
        plt.ylabel('T_B [K]')

    def Tb(self):
        # ####-----Brightness temperature
        plt.figure('brightness')
        for j, b in enumerate(self.data.b):
            Tplt = self.data.Tb[j]
            lt = '-'
            if (len(Tplt) == 1):
                lt = 'o'
            plt.plot(self.data.f, self.data.Tb[j], lt, label=str(b))
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Brightness temperature [K]')
