# -*- mode: python; coding: utf-8 -*-
"""Ray tracing."""
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
import numpy as np
import sys
from . import shape
from . import utils

raypathdir = {'egress': -1, 'ingress': 1, 'tangent': 0}
# Vectors
_X = 0
_Y = 1
_Z = 2
xHat = np.array([1.0, 0.0, 0.0])
yHat = np.array([0.0, 1.0, 0.0])
zHat = np.array([0.0, 0.0, 1.0])


class Ray:
    """Class to hold ray parameters."""

    allowed_parameters = ['ds', 'layer4ds', 'r4ds', 'P4ds', 'doppler', 'tip', 'rotate', 'rNorm']

    def __init__(self):
        """Initialize."""
        for k in self.allowed_parameters:
            setattr(self, k, None)

    def update(self, **kwargs):
        """Update."""
        for k, v, in kwargs.items():
            if k in self.allowed_parameters:
                setattr(self, k, v)
            else:
                raise ValueError('{} not allowed Ray parameter'.format(k))


def computeAspect(Q, f=1.0):
    """Convert the orientation vector [posAng,lat_planetographic] to the rotation angles."""
    tip = -Q[0] * np.pi / 180.0   # 'tip' to north
    # 'rotate' the sub-earth planetographic latitude
    rotate = -np.arctan(np.tan(Q[1] * np.pi / 180.0) * (1.0 - f)**2)
    return tip, rotate   # planetocentric coordinates


def rotate2planet(rotate, tip, b):
    """Separate out to keep use consistent."""
    # out_vec = shape.rotZ(tip,shape.rotX(rotate,b))  # the first way, which is seemingly incorrect
    out_vec = shape.rotX(rotate, shape.rotZ(tip, b))
    return out_vec


def rotate2obs(rotate, tip, b):
    """Opposite to rotate2planet."""
    out_vec = shape.rotZ(tip, shape.rotX(rotate, b))
    return out_vec


def findEdge(atm, b, rNorm, tip, rotate, gtype, printdot=False):
    """Find if at the edge of the planet."""
    tmp = (b[0]**2 + b[1]**2)
    try:
        zQ_Trial = np.arange(np.sqrt(1.0 - tmp) * 1.01, 0.0, -0.005)
    except ValueError:
        return None, None
    pclat_zQ_Trial = []
    r_zQ_Trial = []
    r_pclat_zQ_Trial = []
    hitPlanet = False
    geoid = shape.Shape(gtype)
    for zQ in zQ_Trial:
        if printdot:
            print('.', end='')
            sys.stdout.flush()
        # Q
        b_vec = np.array([b[0], b[1], zQ])
        # --> P
        b_vec = rotate2planet(rotate, tip, b_vec)
        r1 = np.linalg.norm(b_vec) * rNorm
        r_zQ_Trial.append(r1)
        # get planetocentric latitude/longitude
        pclat = utils.r2d(np.arcsin(np.dot(b_vec, yHat) / np.linalg.norm(b_vec)))
        delta_lng = utils.r2d(np.arctan2(np.dot(b_vec, xHat), np.dot(b_vec, zHat)))
        pclat_zQ_Trial.append(pclat)
        r2 = geoid.calcShape(atm, rNorm, pclat, delta_lng)
        r_pclat_zQ_Trial.append(r2)
        if r1 < r2:
            hitPlanet = True
            break
    if not hitPlanet:
        return None, None
    # interpolate to value
    xx = np.flipud(np.array(r_zQ_Trial) - np.array(r_pclat_zQ_Trial))
    yy = np.flipud(np.array(zQ_Trial[0:len(r_zQ_Trial)]))
    try:
        zQ = np.interp(0.0, xx, yy)
        b = np.array([b[0], b[1], zQ])
        edge = rNorm * rotate2planet(rotate, tip, b)
    except ValueError:
        b = None
        edge = None

    del zQ_Trial, pclat_zQ_Trial, r_zQ_Trial, r_pclat_zQ_Trial, geoid, xx, yy
    return edge, b  # same vector but in P and Q coordinate systems


def compute_ds(atm, b, orientation=None, gtype=None, verbose=False):
    """
    Compute the path length through the atmosphere.

    Parameters
    ----------
    b = impact parameter (fractional distance to outer edge at that
    latitude in observer's coordinates) orientation = position angle of
    the planet [0]='tip', [1]='subearth latitude'
    """
    if gtype is None:
        gtype = atm.config.gtype
    if orientation is None:
        orientation = atm.config.orientation
    path = Ray()
    req = atm.property[atm.config.LP['R']]   # radius of layers along equator
    rNorm = req[0]
    nr = atm.property[atm.config.LP['N']]    # refractive index of layers
    if (b[0]**2 + b[1]**2) >= 1.0:
        return path

    mu = np.sqrt(1.0 - b[0]**2 - b[1]**2)

    f = 1.0 - atm.config.Rpol / atm.config.Req
    tip, rotate = computeAspect(orientation, f)
    if verbose:
        print('intersection:  ({:.3f}, {:.3f})    '.format(b[0], b[1]), end='')
        print('aspect:  ({:.4f},  {:.4f})'.format(tip * 180.0 / np.pi, rotate * 180.0 / np.pi))
        print('Finding atmospheric edge', end='')
    edge, b = findEdge(atm, b, rNorm, tip, rotate, gtype)
    if edge is None:
        return path
    r1 = np.linalg.norm(edge)

    # get planetocentric latitude/longitude
    pclat = utils.r2d(np.arcsin(np.dot(edge, yHat) / np.linalg.norm(edge)))
    delta_lng = utils.r2d(np.arctan2(np.dot(edge, xHat), np.dot(edge, zHat)))
    geoid = shape.Shape(gtype)
    r2 = geoid.calcShape(atm, rNorm, pclat, delta_lng)
    if verbose:
        print(' within {:.2f} m'.format(abs(r1 - r2) * 100.0))
        geoid.printShort()

    # initialize - everything below is in planetocentric coordinates, so need to rotate s-vector
    start = np.array([0.0, 0.0, -1.0])
    start = rotate2planet(rotate, tip, start)
    s = [start]
    n = [geoid.n]
    r = [geoid.r]
    t_inc = [np.arccos(-np.dot(s[-1], n[-1]))]            # incident angle_0
    nratio = nr[0] / nr[1]                                # refractive index ratio_0
    t_tran = [np.arcsin(nratio * np.sin(t_inc[-1]))]      # transmitted angle_0

    # return array initialization, use a value to keep indexing numbering for loop
    path.update(ds=[0.0], layer4ds=[0.0], r4ds=[0.0], P4ds=[0.0], doppler=[])

    # loop while in atmosphere
    i = 0                      # loops over the path
    layer = 0                  # keeps track which physical layer you are in
    inAtmosphere = True
    direction = 'ingress'
    while inAtmosphere:
        if verbose:
            print('------------------')
            print('\tstep {}:  layer {} {} '.format(i, layer, direction))
            print('\ts = [{:.4f}, {:.4f}, {:.4f}],  ds = {:.4f}'.
                  format(s[-1][_X], s[-1][_Y], s[-1][_Z], path.ds[-1]))
            geoid.print_a_step()
            print('\tt_inc, tran:  {:.8f} -> {:.8f}'.
                  format(utils.r2d(t_inc[-1]), utils.r2d(t_tran[-1])))
        # update s-vector
        s.append(nratio * s[i] + raypathdir[direction] *
                 (nratio * np.cos(t_inc[i]) * n[i] - np.cos(t_tran[i]) * n[i]))
        rNowMag = geoid.rmag
        rNextMag = geoid.calcShape(atm, req[layer + raypathdir[direction]], pclat, delta_lng)
        rdots = np.dot(r[i], s[i + 1])

        vw = np.interp(pclat, atm.config.vwlat, atm.config.vwdat) / 1000.0
        path.doppler.append(1.0 - (atm.config.omega_m * rNowMag * np.cos(utils.d2r(pclat)) + vw) *
                            np.sin(utils.d2r(delta_lng)) / 3.0E5)

        # get ds, checking for exit, errors and warnings...
        try:
            dsp = -rdots + np.sqrt(rdots**2.0 + rNextMag**2.0 - rNowMag**2.0)
            dsm = -rdots - np.sqrt(rdots**2.0 + rNextMag**2.0 - rNowMag**2.0)
            if direction == 'ingress':
                ds_step = dsm
            elif direction == 'egress':
                ds_step = dsp
        except ValueError:  # tangent layer (probably)
            if direction == 'ingress':
                print('In tangent layer  ', end='')
                direction = 'tangent'
                ds_step = -2.0 * rdots
                if ds_step > rNorm:
                    inAtmosphere = False
                    print('Error:  tangent ds too large')
                    break
                rNextMag = rNowMag
            else:
                inAtmosphere = False
                break
        # What is correct here?  Is this true for egress, or only in error?
        # 8/1/14 I commented out if statement, but only print statement was indented
        if ds_step < 0.0:
            print('Warning:  ds < 0  ({}:  r.s={}, ds={}, [{},{}])'
                  .format(direction, rdots, ds_step, dsp, dsm))
            inAtmosphere = False
            break

        if atm.config.limb == 'sec':  # overwrite ds_step with secant version
            ds_step = abs(rNextMag - rNowMag) / mu

        # append to return arrays
        path.ds.append(ds_step)
        path.layer4ds.append(layer)
        path.r4ds.append(rNowMag)
        path.P4ds.append(atm.gas[atm.config.C['P']][layer])

        # get next step, double-check r value (and compute n, etc)
        rnext = r[i] + path.ds[i + 1] * s[i + 1]
        pclat = utils.r2d(np.arcsin(np.dot(rnext, yHat) / np.linalg.norm(rnext)))
        delta_lng = utils.r2d(np.arctan2(np.dot(rnext, xHat), np.dot(rnext, zHat)))
        r2 = geoid.calcShape(atm, req[layer + raypathdir[direction]], pclat, delta_lng)
        if abs(r2 - rNextMag) > 2.0 and layer != 0:
            print('Warning:  {} != {} ({} km) at layer {}'
                  .format(r2, rNextMag, r2 - rNextMag, layer))
        r.append(rnext)  # which is also geoid.r, or should be
        n.append(geoid.n)

        # get new incident angle
        layer += raypathdir[direction]
        if direction == 'tangent':
            direction = 'egress'
            t_inc.append(np.pi / 2.0)
        else:
            try:
                t_inc.append(np.arccos(-raypathdir[direction] * np.dot(s[i + 1], n[i + 1])))
            except ValueError:
                print('t_inc ValueError |s_(i+1)| = {}, |n_(i+1)| = {} - set to previous'
                      .format(np.linalg.norm(s[i + 1]), np.linalg.norm(n[i + 1])))
                t_inc.append(t_inc[i])
                inAtmosphere = False
                break
        i += 1

        # get refractive index ratio and check for exit
        try:
            nratio = nr[layer] / nr[layer + raypathdir[direction]]
            nratio = 1.0
            try:
                t_tmp = np.arcsin(nratio * np.sin(t_inc[-1]))
            except ValueError:
                t_tmp = nratio * np.pi / 2.0
            t_tran.append(t_tmp)
        except IndexError:
            inAtmosphere = False
        # ## end loop ###

    # Get rid of the first entry, which was just used to make indexing in loop consistent
    del path.ds[0], path.layer4ds[0], path.r4ds[0], path.P4ds[0]

    path.update(tip=tip, rotate=rotate, rNorm=rNorm)

    del s, r, n, geoid, req, nr
    return path
