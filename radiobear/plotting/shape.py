import matplotlib.pyplot as plt
import numpy as np

from .. import shape


# ##############################################################################################################
#                                           SHAPE
# ##############################################################################################################
class plots:
    def __init__(self, shp):
        self.shp = shp

    def show(self):
        plt.show()

    def shapes(self, planet, r, lat=90.0, delta_lng=0.0,
               gtypes=['ellipse', 'sphere', 'reference', 'gravity'], latstep='default'):
        colors = ['k', 'r', 'g', 'b', 'y', 'm', 'c']
        plt.figure('Shapes')
        for i, gtype in enumerate(gtypes):
            self.saveShape = [['oui']]
            self.shp.calcShape(planet, r, pclat=90.0, delta_lng=delta_lng,
                               gtype=gtype, latstep=latstep)
            del self.shp.saveShape[0]  # remove nonsense first term
            self.shp.saveShape = np.array(self.shp.saveShape)
            if gtype != 'gravity':
                self.shp.saveShape = np.flipud(self.shp.saveShape)
            _y = []
            _z = []
            for v in self.shp.saveShape:
                _y.append(v[shape._Y])
                _z.append(v[shape._Z])
            plt.plot(_z, _y, color=colors[self.colorCounter % len(colors)], label=gtype)
            self.colorCounter += 1
            plt.axis('image')
            plt.legend()
