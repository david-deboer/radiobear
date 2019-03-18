import numpy as np
import six
import copy


class DataReturn:
    """
    This holds the data that one may wish to use as output.  Data are stored as numpy arrays
    """

    allowed_parameters = ['f', 'freqUnit', 'b', 'Tb', 'header']

    def __repr__(self):
        s = ''
        for i, b in enumerate(self.b):
            bstr = 'b = {}:'.format(b)
            s += bstr + '\n f = '
            f = ['{:6.1f}'.format(x) for x in self.f]
            s += ' '.join(f) + '  GHz\nTb = '
            T = ['{:6.1f}'.format(x) for x in self.Tb[i]]
            s += ' '.join(T) + '  K\n'
        return s

    def set(self, par, val):
        if par not in self.allowed_parameters:
            print("{} not in valid data return list.".format(par))
            return
        val = copy.copy(val)
        if isinstance(val, list):
            setattr(self, par, np.asarray(val, dtype=np.float32))
        else:
            setattr(self, par, val)

    def show(self):
        for v in self.allowed_parameters:
            if v == 'header':
                self.show_header()
            else:
                print("{}:  {}.".format(v, getattr(self, v)))

    def show_header(self):
        print('Header')
        for k, v in six.iteritems(self.header):
            print("\t{}".format(v.strip()))
