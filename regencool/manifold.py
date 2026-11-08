"""Inlet-manifold flow distribution + per-channel tolerance Monte Carlo.

Header : one inlet at theta = 0 feeding two symmetric branches that run to theta = 180 deg (a ring header). Along a
         branch the static pressure changes by a Bernoulli term for the area change, Darcy friction, and at each
         take-off a branch-entry loss C_b (1/2 rho u_h^2) and a momentum regain k_r (1/2 rho (u_h^2 - u_h'^2)).
         Header area / hydraulic diameter vs. theta are tabulated in the config.
Channel: each channel's response (dP, exit pressure, max Twg, Tout, saturation margin) is tabulated from the 1D
         solver over (feature offset d, flow factor mf) and interpolated bilinearly.
Network: all channels discharge into a common outlet at P_o. For each sample, P_o is found by bisection so that the
         channel flows (each from dP_i(d_i, mf_i) = P_h + off_i - P_o) add up to the total flow.
"""
import math

import numpy as np

from .channel import geometry
from .correlations import f_colebrook
from .solver1d import ChannelClosed


def channel_angles(nch, n_land=0, land_deg=0.0, land0=0.0):
    """Angular positions [deg] of the channel inlets; n_land solid lands of land_deg spaced evenly."""
    if n_land == 0:
        return np.mod(land0 + 360.0 / nch * (np.arange(nch) + 0.5), 360.0)
    pitch = (360.0 - n_land * land_deg) / nch
    counts = [nch // n_land + (1 if b < nch % n_land else 0) for b in range(n_land)]
    out = []
    for b, cnt in enumerate(counts):
        a0 = land0 + b * 360.0 / n_land + land_deg / 2
        out += [a0 + pitch * (k + 0.5) for k in range(cnt)]
    return np.mod(np.array(out), 360.0)


def header_offsets(ang, mdot, rho, mu, hdr, k_r=0.8, C_b=0.5, ks=30e-6, dth=0.25):
    """Static pressure available at each channel entry relative to the inlet static [bar] (uniform withdrawal).
    hdr = dict(theta_deg=[...], area_mm2=[...], dh_mm=[...], R_mean_mm=...)."""
    TH, AR, DHT = (np.asarray(hdr[k], float) for k in ("theta_deg", "area_mm2", "dh_mm"))
    Rm = float(hdr["R_mean_mm"])
    n = len(ang)
    m_ch = mdot / n
    a = np.where(ang > 180, 360 - ang, ang)
    w = np.where(np.isclose(a, 180) | np.isclose(a, 0), 0.5, 1.0)     # channel on a symmetry line: half per branch
    th, Q, p = 0.0, mdot / 2 / rho, 0.0
    order = np.argsort(a)
    done = np.zeros(n, bool)
    off = np.zeros(n)
    prof = []
    A_old = np.interp(0.0, TH, AR) * 1e-6
    while th <= 180.0 + 1e-9:
        A_new = np.interp(th, TH, AR) * 1e-6
        u1, u2 = Q / A_old, Q / A_new
        p -= 0.5 * rho * (u2 * u2 - u1 * u1)
        D = np.interp(th, TH, DHT) * 1e-3
        Re = max(rho * u2 * D / mu, 1.0)
        f = 64 / Re if Re < 2300 else f_colebrook(Re, ks / D)
        p -= f * (Rm * 1e-3 * math.radians(dth)) / D * 0.5 * rho * u2 * u2
        for i in order:
            if not done[i] and a[i] <= th + 1e-9:
                uh = Q / A_new
                off[i] = p - C_b * 0.5 * rho * uh * uh
                Qn = max(Q - w[i] * m_ch / rho, 0.0)
                un = Qn / A_new
                p += k_r * 0.5 * rho * (uh * uh - un * un)
                Q = Qn
                done[i] = True
        prof.append((th, p / 1e5, Q / A_new))
        A_old = A_new
        th += dth
    return off / 1e5, np.array(prof)


def tables(case, d_grid, m_grid, C_hg=None, log=print):
    """Channel response tables [len(d_grid), len(m_grid)] (flow factor relative to nominal per-channel flow)."""
    keys = ("dP", "dPex", "Twg", "Tout", "mrg")
    T = {k: np.full((len(d_grid), len(m_grid)), np.nan) for k in keys}
    fl = case.fluid
    Pin_b = case.Pin / 1e5
    for i, d in enumerate(d_grid):
        try:
            D = case.design(d)
            G = geometry(D, case.contour, case.dx)
        except (ChannelClosed, ValueError):
            continue
        for j, mf in enumerate(m_grid):
            try:
                r = case.solve(D, mdot=case.mdot * mf, C_hg=C_hg, G=G)
            except (ChannelClosed, ValueError):
                continue
            mrg = 99.0
            if None not in (fl.Tc, fl.Pc, fl.omega):
                sub = r["Tb"] < fl.Tc
                if sub.any():
                    mrg = float(np.min(r["P"][sub] / 1e5 - fl.psat_bar(r["Tb"][sub])))
            T["dP"][i, j] = r["dP_tot"]
            T["dPex"][i, j] = Pin_b - r["Pmin"]
            T["Twg"][i, j] = r["Twg"].max() - 273.15
            T["Tout"][i, j] = r["Tout"] - 273.15
            T["mrg"][i, j] = mrg
        log("  table d = %+.3f mm done" % d)
    for k in keys:                      # fill failures (closed channel / tiny flow) conservatively
        a = T[k]
        for i in range(a.shape[0]):
            bad = np.isnan(a[i])
            if bad.all():
                a[i] = 1e3 if k in ("dP", "dPex", "Twg") else -1e3
            elif bad.any():
                a[i, bad] = np.interp(np.nonzero(bad)[0], np.nonzero(~bad)[0], a[i, ~bad])
    T["d"], T["mf"] = np.asarray(d_grid, float), np.asarray(m_grid, float)
    return T


def interp2(T, key, d, mf):
    DG, MG, tab = T["d"], T["mf"], T[key]
    x = (np.clip(d, DG[0], DG[-1]) - DG[0]) / (DG[1] - DG[0])
    y = (np.clip(mf, MG[0], MG[-1]) - MG[0]) / (MG[1] - MG[0])
    i = np.minimum(np.floor(x).astype(int), len(DG) - 2)
    j = np.minimum(np.floor(y).astype(int), len(MG) - 2)
    fx, fy = x - i, y - j
    return (tab[i, j] * (1 - fx) * (1 - fy) + tab[i + 1, j] * fx * (1 - fy) + tab[i, j + 1] * (1 - fx) * fy
            + tab[i + 1, j + 1] * fx * fy)


def network(T, d, off, P_h, P_table=None):
    """d [S, n] channel offsets [mm], off [n] header offsets [bar], P_h inlet static [bar].
    Returns per-channel flow factors mf (sum = n per sample) and the interpolated responses."""
    P_table = P_h if P_table is None else P_table
    S, n = d.shape
    MG = T["mf"]

    def m_of(dPav):
        lo, hi = np.full((S, n), MG[0]), np.full((S, n), MG[-1])
        for _ in range(32):
            mid = 0.5 * (lo + hi)
            g = interp2(T, "dP", d, mid) - dPav
            hi = np.where(g > 0, mid, hi)
            lo = np.where(g > 0, lo, mid)
        return 0.5 * (lo + hi)

    dmax = float(np.nanmax(np.where(T["dP"] < 1e3, T["dP"], np.nan)))
    plo, phi = np.full(S, P_h - 2 * dmax - 1), np.full(S, P_h)
    for _ in range(36):
        po = 0.5 * (plo + phi)
        tot = m_of(P_h + off[None, :] - po[:, None]).sum(1)
        phi = np.where(tot < n, po, phi)
        plo = np.where(tot < n, plo, po)
    po = 0.5 * (plo + phi)
    mf = m_of(P_h + off[None, :] - po[:, None])
    Pex = P_h + off[None, :] - interp2(T, "dPex", d, mf)
    return dict(mf=mf, Po=po, Pex=Pex, Twg=interp2(T, "Twg", d, mf), Tout=interp2(T, "Tout", d, mf),
                mrg=interp2(T, "mrg", d, mf) - (P_table - (P_h + off[None, :])))


def monte_carlo(T, off, P_h, sig, mu=0.0, S=2000, seed=1):
    n = len(off)
    r = np.random.default_rng(seed)
    d = mu + sig * r.standard_normal((S, n)) if sig > 0 else np.full((S, n), mu)
    R = network(T, d, off, P_h)

    def q(a, p):
        return float(np.percentile(a, p))

    mfmin, Tmax = R["mf"].min(1), R["Twg"].max(1)
    return dict(mu=mu, sigma=sig, samples=S, mf_min_p50=q(mfmin, 50), mf_min_p01=q(mfmin, 1),
                Twg_max_p50=q(Tmax, 50), Twg_max_p99=q(Tmax, 99), Pex_min_p01=q(R["Pex"].min(1), 1),
                psat_margin_p01=q(R["mrg"].min(1), 1), Po_p50=q(R["Po"], 50), Tout_max_p99=q(R["Tout"].max(1), 99),
                _Tmax=Tmax, _mfmin=mfmin)
