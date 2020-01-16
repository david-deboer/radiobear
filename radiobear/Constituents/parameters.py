from argparse import Namespace


def setpar(kwargs):
    par = Namespace(units='dBperkm', path='./', verbose=False)
    for p, v in kwargs.items():
        setattr(par, p, v)
    return par
