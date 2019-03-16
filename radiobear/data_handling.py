class DataReturn:
    f = []
    b = []
    Tb = []
    header = {}

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
