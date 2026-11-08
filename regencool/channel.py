"""Cooling-channel design and its geometry along the contour.

All lengths in mm. t(x): hot-wall thickness and s(x): channel size are measured normal to the contour.
s is the height for 'rect', the diameter for 'circ' and the total height (bottom to apex) for 'tear'.
'tear' = circle of radius R = s/(1+sqrt 2) with a 45-deg roof (printable without support in LPBF).
"""
import copy
import math

import numpy as np
from scipy.interpolate import PchipInterpolator

TEAR_A, TEAR_P, TEAR_H = 0.75 * math.pi + 1.0, 1.5 * math.pi + 2.0, 1.0 + math.sqrt(2.0)   # x R^2, x R, x R
SHAPES = ("rect", "circ", "tear")


def knots_fun(knots, log=False):
    """Piecewise-cubic monotone (PCHIP) function from [[x, value], ...], held constant outside the range."""
    k = np.asarray(knots, float).reshape(-1, 2)
    k = k[np.argsort(k[:, 0])]
    if len(k) == 1:
        v = float(k[0, 1])
        return lambda x: np.full_like(np.asarray(x, float), v)
    y = np.log(k[:, 1]) if log else k[:, 1]
    P = PchipInterpolator(k[:, 0], y)
    lo, hi = k[0, 0], k[-1, 0]
    if log:
        return lambda x: np.exp(P(np.clip(np.asarray(x, float), lo, hi)))
    return lambda x: P(np.clip(np.asarray(x, float), lo, hi))


class ChannelDesign:
    """nch channels, optional n_land solid lands of land_deg each (e.g. instrumentation bosses).

    beta_fun(x): helix angle [deg] of the channel path (0 = axial). Rect channels can be split by a rib of
    width rib_s downstream of x >= xs, and are split automatically while the sub-channel width exceeds w_split.
    """

    def __init__(self, t_fun, s_fun, nch, shape="tear", rib=1.0, land_deg=0.0, n_land=0, beta_fun=None,
                 xs=None, rib_s=None, w_split=None, name=""):
        if shape not in SHAPES:
            raise ValueError("shape must be one of %s" % (SHAPES,))
        self.t, self.s, self.nch, self.shape, self.rib = t_fun, s_fun, int(nch), shape, rib
        self.land, self.NL = land_deg, n_land
        self.beta_fun, self.xs, self.rib_s, self.w_split, self.name = beta_fun, xs, rib_s, w_split, name
        self.Neff = self.nch * 360.0 / (360.0 - n_land * land_deg)   # equivalent pitch count over 360 deg
        self.dh = 0.0    # feature-size offset (tolerance studies)
        self.dw = 0.0    # rib/web offset (rect)

    def copy(self):
        return copy.copy(self)

    def offset(self, d):
        """Manufacturing offset: channel features +d, rect ribs -d, hot wall -d/2 (material removed)."""
        D = self.copy()
        t0 = self.t
        D.t = lambda x, t0=t0: t0(x) - d / 2
        D.dh = self.dh + d
        D.dw = self.dw + (d if self.shape == "rect" else 0.0)
        return D

    def clamp_web(self, contour, web_min, dx=0.25):
        """circ/tear: limit s(x) so that the web between neighbouring channels stays >= web_min."""
        if self.shape == "rect" or not web_min:
            return self
        G = geometry(self, contour, dx)
        smax = G["pitch"] - web_min
        if self.shape == "circ":
            sc = np.minimum(G["s"], np.maximum(smax, 0.3))
        else:
            sc = np.minimum(G["s"], np.maximum(smax / 2 * TEAR_H, 0.3))
        xg = G["x"]
        D = self.copy()
        D.s = lambda x, xg=xg, sc=sc: np.interp(x, xg, sc)
        return D


def geometry(D, contour, dx=0.25):
    """Arrays along x: wall radius r, slope, channel mid radius rm, helix angle beta [deg], normal pitch,
    path length per unit x (ds), normal curvature kn and geodesic curvature kg [1/mm]."""
    xg = np.arange(contour.x0, contour.x_end + 1e-9, dx)
    r = contour.R(xg)
    rp = contour.dR(xg)
    ca = 1 / np.sqrt(1 + rp ** 2)
    t = np.asarray(D.t(xg), float) * np.ones_like(xg)
    s = np.asarray(D.s(xg), float) * np.ones_like(xg)
    rm = r + (t + s / 2) * ca
    rmp = np.gradient(rm, xg)
    rmpp = np.gradient(rmp, xg)
    smer = np.sqrt(1 + rmp ** 2)
    beta_deg = np.zeros_like(xg) if D.beta_fun is None else np.asarray(D.beta_fun(xg), float) * np.ones_like(xg)
    thp = np.tan(np.radians(beta_deg)) * smer / rm          # d(theta)/dx [rad/mm]
    ds = np.sqrt(smer ** 2 + (rm * thp) ** 2)
    cosb, sinb = smer / ds, rm * thp / ds
    pitch = 2 * math.pi * rm / D.Neff * cosb
    k1 = rmpp / (1 + rmp ** 2) ** 1.5
    k2 = -1 / (rm * np.sqrt(1 + rmp ** 2))
    kn = k1 * cosb ** 2 + k2 * sinb ** 2
    beta = np.arctan2(sinb, cosb)
    kg = np.gradient(beta, xg) / ds
    return dict(x=xg, r=r, rp=rp, ca=ca, t=t, s=s, rm=rm, ds=ds, cosb=cosb, pitch=pitch, kn=kn, kg=kg,
                beta=np.degrees(beta))


def section(D, G, j):
    """Cross-section at grid index j.

    ns: sub-channels per pitch, we/hh: width/height, A: area, Dh: hydraulic diameter, Pd: wetted perimeter
    facing the hot wall directly, Pf: fin (rib side) perimeter, Hf: fin height, tf: fin thickness,
    teff: effective conduction thickness, a: aspect ratio, web: minimum solid between channels.
    """
    pitch, s, t = G["pitch"][j], G["s"][j] + D.dh, G["t"][j]
    if D.shape == "rect":
        w = pitch - (D.rib - D.dw)
        ns = 1
        if D.xs is not None and G["x"][j] >= D.xs:
            ns = 2
        rs = (D.rib_s if D.rib_s is not None else D.rib) - D.dw
        if D.w_split:
            while (w - (ns - 1) * rs) / ns > D.w_split + D.dh:
                ns += 1
        we = (w - (ns - 1) * rs) / ns
        hh = s
        A = we * hh
        return dict(ns=ns, we=we, hh=hh, A=A, Dh=2 * we * hh / (we + hh), Pd=we, Pf=2 * hh, Hf=hh,
                    tf=D.rib - D.dw, teff=t, a=hh / we, web=D.rib - D.dw, rs=rs, P=2 * (we + hh))
    if D.shape == "circ":
        d = s
        A = math.pi * d * d / 4
        return dict(ns=1, we=d, hh=d, A=A, Dh=d, Pd=math.pi * d / 2, Pf=math.pi * d / 2, Hf=d / 2,
                    tf=max(pitch - math.pi * d / 4, 0.2), teff=t + 0.107 * d, a=1.0, web=pitch - d, rs=0.0,
                    P=math.pi * d)
    R = s / TEAR_H
    A = TEAR_A * R * R
    P = TEAR_P * R
    Aup = A - 0.5 * math.pi * R * R
    Hf = math.sqrt(2) * R
    return dict(ns=1, we=2 * R, hh=s, A=A, Dh=4 * A / P, Pd=math.pi * R, Pf=P - math.pi * R, Hf=Hf,
                tf=max(pitch - Aup / Hf, 0.2), teff=t + 0.107 * 2 * R, a=1.0, web=pitch - 2 * R, rs=0.0, P=P)
