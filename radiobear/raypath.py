# -*- mode: python; coding: utf-8 -*-
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
from __future__ import absolute_import, division, print_function
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
    def __init__(self, ds=None, layer4ds=None, r4ds=None, P4ds=None, doppler=None,
                 tip=None, rotate=None, rNorm=None):
        self.ds = ds
        self.layer4ds = layer4ds
        self.r4ds = r4ds
        self.P4ds = P4ds
        self.rNorm = rNorm
        self.tip = tip
        self.rotate = rotate
        self.doppler = doppler

    def update(self, ds=None, layer4ds=None, r4ds=None, P4ds=None, doppler=None, tip=None,
               rotate=None, rNorm=None):
        if ds is not None:
            self.ds = ds
        if layer4ds is not None:
            self.layer4ds = layer4ds
        if r4ds is not None:
            self.r4ds = r4ds
        if P4ds is not None:
            self.P4ds = P4ds
        if rNorm is not None:
            self.rNorm = rNorm
        if tip is not None:
            self.tip = tip
        if rotate is not None:
            self.rotate = rotate
        if doppler is not None:
            self.doppler = doppler


def computeAspect(Q, f=1.0):
    """Convert the orientation vector [posAng,lat_planetographic] to the rotation angles"""
    tip = -Q[0] * np.pi / 180.0   # 'tip' to north
    # 'rotate' the sub-earth planetographic latitude
    rotate = -np.arctan(np.tan(Q[1] * np.pi / 180.0) * (1.0 - f)**2)
    return tip, rotate   # planetocentric coordinates


def rotate2planet(rotate, tip, b):
    """separate out to keep use consistent!"""
    # out_vec = shape.rotZ(tip,shape.rotX(rotate,b))  # the first way, which is seemingly incorrect
    out_vec = shape.rotX(rotate, shape.rotZ(tip, b))
    return out_vec


def rotate2obs(rotate, tip, b):
    """This should be opposite to rotate2planet..."""
    out_vec = shape.rotZ(tip, shape.rotX(rotate, b))
    return out_vec


def findEdge(atm, b, rNorm, tip, rotate, gtype, printdot=False):
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
    """Computes the path length through the atmosphere given:
            b = impact parameter (fractional distance to outer edge at that
            latitude in observer's coordinates) orientation = position angle of
            the planet [0]='tip', [1]='subearth latitude' """
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
    t_tran = [np.arcsin(nratio * np.sin(t_inc[-1]))]    # transmitted angle_0

    # return array initialization, use a value to keep indexing numbering for loop
    ds = [0.0]
    layer4ds = [0.0]
    r4ds = [0.0]  # not needed, but for information
    P4ds = [0.0]  # not needed, but for information
    doppler = []

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
                  format(s[-1][_X], s[-1][_Y], s[-1][_Z], ds[-1]))
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
        dopp = (1.0 - (atm.config.omega_m * rNowMag * np.cos(utils.d2r(pclat)) + vw) *
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
            print('Error:  ds < 0  ({}:  r.s={}, ds={}, [{},{}])'
                  .format(direction, rdots, ds_step, dsp, dsm))
            inAtmosphere = False
            break

        if atm.config.limb == 'sec':  # overwrite ds_step with secant version
            ds_step = abs(rNextMag - rNowMag) / mu

        # append to return arrays
        ds.append(ds_step)
        layer4ds.append(layer)
        r4ds.append(rNowMag)
        P4ds.append(atm.gas[atm.config.C['P']][layer])
        doppler.append(dopp)

        # get next step, double-check r value (and compute n, etc)
        rnext = r[i] + ds[i + 1] * s[i + 1]
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
    del ds[0], layer4ds[0], r4ds[0], P4ds[0]
    dsmu = []
    for i in range(min(len(ds), len(r4ds)) - 1):
        if abs(r4ds[i] - r4ds[i + 1]) < 1E-6:
            dsmuappend = 0.001
        else:
            dsmuappend = ds[i] / (r4ds[i] - r4ds[i + 1])
        dsmu.append(dsmuappend)
    path.update(ds=ds, layer4ds=layer4ds, r4ds=r4ds, P4ds=P4ds, doppler=doppler,
                tip=tip, rotate=rotate, rNorm=rNorm)

    del s, r, n, ds, layer4ds, r4ds, P4ds, geoid, req, nr
    return path


# ##Test functions
def refractTest(layers):
    """Returns fake refractive index at layers"""
    n = [1.0]
    n = []
    for r in layers:
        v = 1.0 + (3000.0 / r)**2
        n.append(v)
    return n


def computeEdges(z):  # not used anymore
    edge = [(3.0 * z[0] - z[1]) / 2.0]
    for i in range(len(z) - 1):
        edge.append((z[i] + z[i + 1]) / 2.0)
    edge.append((3.0 * z[-1] - z[-2]) / 2.0)
    return edge


def layersTest(rmin=100.0, rmax=20000.0, nlyr=100):
    """Returns layer"""
    dr = (rmax - rmin) / nlyr
    mid = []
    for i in range(nlyr):
        r = rmax - i * dr
        mid.append(r)
    return mid


def testPath(b=0.5, rmin=12000.0, rmax=20000.0, nlyr=50, verbose=False):
    # make layers
    mid = layersTest(rmin=rmin, rmax=rmax, nlyr=nlyr)
    n = refractTest(mid)
    ds = compute_ds(b, mid, n, verbose=verbose)
    return mid, ds
