from argparse import Namespace
import six


def setpar(kwargs):
    par = Namespace(units='dBperkm', path='./', verbose=False)
    for p, v in six.iteritems(kwargs):
        setattr(par, p, v)
    return par
