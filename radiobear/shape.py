from __future__ import absolute_import, division, print_function
import numpy as np
import six
import scipy.special as scisp
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from . import utils as u
from . import atmosphere
_X = 0
_Y = 1
_Z = 2


class Shape:
    """Calculates radius and normal to planet.  Types are:
            gravity, reference, sphere/circle, ellipse"""
    def __init__(self, gtype='ellipse'):
        self.gtype = gtype
        self.gamma = 0.0
        self.pclat = 0.0
        self.pglat = 0.0
        self.delta_lng = 0.0
        self.omega = 0.0
        self.r = np.zeros(3)
        self.rmag = 0.0
        self.n = np.zeros(3)
        self.t = np.zeros(3)
        self.g = []
        self.gmag = 0.0
        self.g_static = 0.0
        self.default_gravcalc_latstep = 0.01  # it needs to be this small to be accurate, but then the gravity calc is too slow
        self.referenceGeoid = False
        self.referenceRadius = 0.0
        self.calcCounter = 0
        self.colorCounter = 0
        self.saveShape = False  # This is the one concession to get a gravity profile without looping many times...so fake it for ellipse

    def __str__(self):
        s = '-------------------------Geoid---------------------------\n'
        s += 'type:  %s, gamma = %.2f deg, lat = %.2f deg, lng = %.2f\n' % (self.gtype, u.r2d(self.gamma), u.r2d(self.pclat), u.r2d(self.delta_lng))
        s += 'r = [%.2f,  %.2f,  %.2f] km, |r| = %.2f km\n' % (self.r[_X], self.r[_Y], self.r[_Z], self.rmag)
        s += 'n = [%.4f,  %.4f,  %.4f]\n' % (self.n[_X], self.n[_Y], self.n[_Z])
        s += 't = [%.4f,  %.4f,  %.4f]\n' % (self.t[_X], self.t[_Y], self.t[_Z])
        s += 'g = [%.6f,  %.6f], |g| = %.6f\n' % (self.g[0], self.g[1], self.gmag)
        s += '---------------------------------------------------------\n'
        return s

    def printShort(self):
        s = '%s:  gamma = %.2f, lat = %.2f, lng = %.2f, |r|=%.2f km\n' % (self.gtype, u.r2d(self.gamma), u.r2d(self.pclat), u.r2d(self.delta_lng), self.rmag)
        print(s)
        return s

    def print_a_step(self):
        s = '\ttype:  %s, gamma = %.2f deg, lat = %.2f deg, lng = %.2f  |g| = %.2f\n' % (self.gtype, u.r2d(self.gamma), u.r2d(self.pclat), u.r2d(self.delta_lng), self.gmag * 1000.0)
        s += '\tr = [%.2f,  %.2f,  %.2f] km, |r| = %.2f km\n' % (self.r[_X], self.r[_Y], self.r[_Z], self.rmag)
        s += '\tn = [%.4f,  %.4f,  %.4f]' % (self.n[_X], self.n[_Y], self.n[_Z])
        print(s)
        return s

    # ##------------------------------------General handling function----------------------------------------
    def calcShape(self, planet, r, pclat=90.0, delta_lng=0.0, gtype=None, latstep='default'):
        # set/reset gtype and latstep (ellipse doesn't need latstep)
        if gtype is None:
            gtype = self.gtype
        self.gtype = gtype
        if isinstance(latstep, six.string_types):
            self.gravcalc_latstep = self.default_gravcalc_latstep
        else:
            self.gravcalc_latstep = float(latstep)

        if gtype == 'gravity':
            r = self.__calcGeoid(planet, r, pclat=pclat, delta_lng=delta_lng)
        elif gtype == 'reference':
            r = self.__calcFromReference(planet, r, pclat=pclat, delta_lng=delta_lng)
        elif gtype == 'ellipse' or gtype == 'circle' or gtype == 'sphere':
            r = self.__calcEllipse(planet, r, pclat=pclat, delta_lng=delta_lng)
        else:
            print(gtype + ' not a valid planet shape')
            r = None
        return r

    # ##---------------------Below are the specific shape handlers----------------------
    # ##'reference' - fits a geoid at the reference pressure and scales from there
    def __calcFromReference(self, planet, r, pclat=90.0, delta_lng=0.0):
        if isinstance(self.referenceGeoid, bool):
            print('Calculating reference geoid at P_ref={}'.format(planet.config.p_ref))
            self.referenceRadius = planet.config.Req
            self.__calcGeoid(planet, planet.config.Req, 90, 0.0)
            self.referenceGeoid = np.array(self.referenceGeoid)
        rlat = np.interp(pclat, self.referenceGeoid[:, 0], self.referenceGeoid[:, 1]) * (r / self.referenceRadius)
        gamma = np.interp(pclat, self.referenceGeoid[:, 0], self.referenceGeoid[:, 2])

        lat = u.d2r(pclat)
        lng = u.d2r(delta_lng)
        norm = np.array([0.0, np.sin(lat + gamma), np.cos(lat + gamma)])
        norm = rotY(lng, norm)
        tang = np.array([0.0, np.cos(lat + gamma), -np.sin(lat + gamma)])
        tang = rotY(lng, tang)
        r_vec = np.array([0.0, rlat * np.sin(lat), rlat * np.cos(lat)])
        r_vec = rotY(lng, r_vec)
        self.gamma = gamma
        self.pclat = lat
        self.pglat = lat + gamma
        self.delta_lng = lng
        self.r = r_vec
        self.rmag = np.linalg.norm(r_vec)
        self.n = norm
        self.t = tang
        self.g = [0.0, 0.0]
        self.gmag = 0.0

        if self.saveShape:  # This is here to fake saveShape to act like __calcGeoid
            plot_latstep = 0.1  # can't be as small as gravcalc since too many recursions
            self.saveShape.append(list(self.r))
            pclat -= plot_latstep
            if pclat > 0.0:
                self.__calcFromReference(planet, r, pclat, delta_lng)
        return self.rmag

    # ##'gravity' - does the full thing, but is very time-consuming
    def __calcGeoid(self, planet, r, pclat, delta_lng):
        """Starts at equatorial radius and moves north or south to compute geoid at pclat"""
        self.calcCounter += 1
        print('Geoid calc counter: ', self.calcCounter)

        if pclat == 0.0:
            nsp = 1.0
        else:
            nsp = np.sign(pclat)
        latstep = nsp * self.gravcalc_latstep
        pclatSteps = np.arange(0.0, pclat + latstep, latstep)

        if self.gtype == 'reference' and type(self.referenceGeoid) == bool:
            self.referenceGeoid = []

        GM = np.interp(r, planet.layerProperty[planet.config.LP['R']], planet.layerProperty[planet.config.LP['GM']])
        for latv in pclatSteps:
            radlatv = u.d2r(latv)
            vw = np.interp(latv, planet.config.vwlat, planet.config.vwdat) / 1000.0
            self.omega = planet.config.omega_m + vw / (r * np.cos(radlatv))
            self.__gravity__(latv, delta_lng, r, GM, self.omega, planet.config.Jn, planet.config.RJ)
            dlat = u.d2r(latstep)
            dr_vec = r * dlat * self.t
            r_vec = self.r + dr_vec
            r = np.linalg.norm(r_vec)
            if self.gtype == 'reference':
                self.referenceGeoid.append([latv, r, self.gamma])
            if self.saveShape:
                self.saveShape.append(list(self.r))
        del pclatSteps
        return self.rmag

    def __gravity__(self, pclat, delta_lng, r, GM, omega, Jn, RJ):
        self.g_static = GM / r**2
        lat = u.d2r(pclat)
        lng = u.d2r(delta_lng)
        if lat == 0.0:
            nsl = 1.0
        else:
            nsl = np.sign(lat)
        dphi = nsl * 0.00001
        Sr = 0.0
        Sp = 0.0
        sp = np.sin(lat)
        sp1 = np.sin(lat + dphi)
        sp0 = np.sin(lat - dphi)
        for i in range(len(Jn)):
            # g(r)
            Sr += (i + 1.0) * Jn[i] * pow(RJ / r, i) * scisp.legendre(i)(sp)
            # g(phi)
            dP = (scisp.legendre(i)(sp1) - scisp.legendre(i)(sp)) / dphi
            dP += (scisp.legendre(i)(sp) - scisp.legendre(i)(sp0)) / dphi
            dP *= 0.5
            Sp += Jn[i] * pow(RJ / r, i) * dP
        # g(r)
        gr = self.g_static * (1.0 - Sr) - (2.0 / 3.0) * (omega**2.0) * r * (1.0 - scisp.legendre(2)(sp))
        # g(phi)
        dP = (3.0 * sp * np.sqrt(1.0 - sp**2))
        gp = (1.0 / 3.0) * (omega**2.0) * r * dP + self.g_static * Sp
        gt = np.sqrt(gr**2 + gp**2)
        # geoid
        gamma = np.arctan2(gp, gr)
        norm = np.array([0.0, np.sin(lat + gamma), np.cos(lat + gamma)])
        norm = rotY(lng, norm)
        tang = np.array([0.0, np.cos(lat + gamma), -np.sin(lat + gamma)])
        tang = rotY(lng, tang)
        r_vec = np.array([0.0, r * np.sin(lat), r * np.cos(lat)])
        r_vec = rotY(lng, r_vec)
        self.gamma = gamma
        self.pclat = lat
        self.pglat = lat + gamma
        self.delta_lng = lng
        self.r = r_vec
        self.rmag = np.linalg.norm(r_vec)
        self.n = norm
        self.t = tang
        self.g = [gr, gp]
        self.gmag = gt

    # ## 'ellipse' or 'circle' or 'sphere' - simply does the equivalent ellipse or circle
    def __calcEllipse(self, planet, r, pclat, delta_lng):
        a = r
        if self.gtype == 'ellipse':
            b = (planet.config.Rpol / planet.config.Req) * r
        else:
            b = r
        lat = u.d2r(pclat)
        lng = u.d2r(delta_lng)
        if lat == 0.0:
            nsl = 1.0
            lat = 1.0E-6
        else:
            nsl = np.sign(lat)

        norm = np.array([0.0, a * np.sin(lat), b * np.cos(lat)])
        norm = norm / np.linalg.norm(norm)
        norm = rotY(lng, norm)
        tang = np.array([0.0, -b * np.cos(lat), a * np.sin(lat)])
        tang = tang / np.linalg.norm(tang)
        tang = rotY(lng, tang)
        r_vec = np.array([0.0, b * np.sin(lat), a * np.cos(lat)])
        r_vec = rotY(lng, r_vec)
        self.rmag = np.linalg.norm(r_vec)
        GM = np.interp(r, planet.layerProperty[planet.config.LP['R']], planet.layerProperty[planet.config.LP['GM']])
        self.g_static = GM / self.rmag**2
        try:
            self.gamma = nsl * np.arccos(np.dot(r_vec, norm) / self.rmag)  # don't need to worry about direction here, since norm etc defined
        except ValueError:
            arg = np.dot(r_vec, norm) / self.rmag
            print('gamma warning (%s):  [{:.2f},{:.2f},{:.2f},{:.2f}]'.format(self.gtype, arg, a, b, pclat), end='')
            print('...but proceeding anyway by setting gamma=0.0')
            self.gamma = 0.0
        self.pclat = lat
        self.pglat = lat + self.gamma
        self.delta_lng = lng
        self.r = r_vec
        self.n = norm
        self.t = tang
        self.g = [self.g_static, self.g_static]
        self.gmag = self.g_static

        if self.saveShape:  # This is here to fake saveShape to act like __calcGeoid
            plot_latstep = 0.1  # can't be as small as gravcalc since too many recursions
            self.saveShape.append(list(self.r))
            pclat -= plot_latstep
            if pclat > 0.0:
                self.__calcEllipse(planet, r, pclat, delta_lng)

        return self.rmag

    def plotShapes(self, planet, r, lat=90.0, delta_lng=0.0, gtypes=['ellipse', 'sphere', 'reference', 'gravity'], latstep='default'):
        colors = ['k', 'r', 'g', 'b', 'y', 'm', 'c']
        plt.figure('Shapes')
        for i, gtype in enumerate(gtypes):
            self.saveShape = [['oui']]
            rmag = self.calcShape(planet, r, pclat=90.0, delta_lng=delta_lng, gtype=gtype, latstep=latstep)
            del self.saveShape[0]  # remove nonsense first term
            self.saveShape = np.array(self.saveShape)
            if gtype != 'gravity':
                self.saveShape = np.flipud(self.saveShape)
            _y = []
            _z = []
            for v in self.saveShape:
                _y.append(v[_Y])
                _z.append(v[_Z])
            plt.plot(_z, _y, color=colors[self.colorCounter % len(colors)], label=gtype)
            self.colorCounter += 1
            plt.axis('image')
            plt.legend()


def rotX(x, V):
    Rx = np.array([[1.0, 0.0, 0.0],
                   [0.0, np.cos(x), -np.sin(x)],
                   [0.0, np.sin(x), np.cos(x)]])
    return np.dot(Rx, V)


def rotY(y, V):
    Ry = np.array([[np.cos(y), 0.0, np.sin(y)],
                   [0.0, 1.0, 0.0],
                   [-np.sin(y), 0.0, np.cos(y)]])
    return np.dot(Ry, V)


def rotZ(z, V):
    Rz = np.array([[np.cos(z), -np.sin(z), 0.0],
                   [np.sin(z), np.cos(z), 0.0],
                   [0.0, 0.0, 1.0]])
    return np.dot(Rz, V)
