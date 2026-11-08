"""2D thermo-elastic check of one channel sector at an axial station, driven by a 1D solution.

Boundary data from the 1D march at that x: hg, Taw (gas side), hc, Tb, P (coolant side), p_gas.
hc is rescaled so that hc_eff x (meshed wetted perimeter) = hc x (true perimeter / cos beta): the 2D mesh resolves
fin conduction itself, so no fin efficiency is applied. Axial end load per sector:
    Fz = p_gas pi r_in^2 / N_sector + P_cool A_channel / cos(beta)
Rect channels split into n_sub sub-channels are modelled as one sub-channel sector (2 pi / (Neff n_sub)).
Best used where the wall slope is small (cylinder, throat); the section is a transverse cut at fixed x.
"""
import math
import warnings

import numpy as np

from . import fe2d
from .channel import section


def station(x, D, res, mat, jacket=2.0, h=0.015, t=None, keep_field=False):
    """x [mm], D: ChannelDesign used for res (a solver1d result), mat: fe2d.Material, jacket: closeout thickness
    above the channel [mm], h: element size near the hot wall [mm], t: optional hot-wall override [mm]."""
    G = res["G"]
    j = int(np.argmin(np.abs(G["x"] - x)))
    x = float(G["x"][j])
    if abs(G["rp"][j]) > 0.25:
        warnings.warn("station x=%.1f: wall slope %.2f, the transverse-section model is approximate" % (x, G["rp"][j]))
    S = section(D, G, j)
    r_in = float(G["r"][j])
    s = float(G["s"][j] + D.dh)
    beta = float(G["beta"][j])
    t = float(G["t"][j]) if t is None else float(t)
    ns = int(S["ns"])
    dth = 2 * math.pi / D.Neff / ns
    hg, Taw, hc, Tb = (float(res[k][j]) for k in ("hg", "Taw", "hc", "Tb"))
    Pc, pg = float(res["P"][j]) / 1e6, float(res["pg"][j]) / 1e6           # MPa
    r_out = r_in + t + s + jacket
    M = fe2d.mesh(r_in, r_out, dth, t, D.shape, s, beta, we=S["we"], h=h)
    Ed = fe2d.edges(M)
    cb = math.cos(math.radians(beta))
    hc_eff = hc * (S["P"] / cb) / Ed["cool"][:, 4].sum()
    T, _ = fe2d.thermal(M, hg, Taw, hc_eff, Tb, mat.k)
    Fz = pg * math.pi * r_in ** 2 / (D.Neff * ns) + Pc * S["A"] / cb
    st, geo = fe2d.structural(M, T, Pc, pg, Fz, mat)
    gasn = np.unique(Ed["gas"][:, :2].astype(int))
    Tc = T[M["conn"]].mean(1)
    cn = np.unique(Ed["cool"][:, :2].astype(int))
    adj = np.isin(M["conn"], cn).any(1)                       # elements touching the channel wall
    inn = M["I"] == 0                                         # hot-wall surface row
    mid = inn & (M["J"] == M["nt"] // 2)                      # under the channel centre
    web = inn & ((M["J"] == 0) | (M["J"] == M["nt"] - 1))     # under the web (sector edge)
    jack = geo["rc"] > r_in + M["apex"] + 0.05
    hot = geo["rc"] < r_in + t
    out = dict(x=x, t=t, s=s, beta=beta, n_sub=ns, r_in=r_in, r_out=r_out, jacket=jacket, hg=hg, Taw=Taw - 273.15,
               hc=hc, hc_eff=hc_eff, Tb=Tb - 273.15, P_cool=Pc * 10, p_gas=pg * 10, n_elem=len(M["conn"]),
               Twg_max=float(np.nanmax(T[gasn])) - 273.15, Twg_min=float(np.nanmin(T[gasn])) - 273.15,
               Twg_1d=float(res["Twg"][j]) - 273.15,
               T_jacket=float(np.average(Tc[jack], weights=geo["area"][jack])) - 273.15 if jack.any() else float("nan"),
               T_hotwall_mean=float(np.average(Tc[hot], weights=geo["area"][hot])) - 273.15)
    for case, d in st.items():
        k = int(np.argmax(d["vm"]))
        k2 = int(np.argmax(np.where(adj, -1, d["vm"])))
        out[case] = dict(vm_max=float(d["vm"][k]), vm_max_dr=float(geo["rc"][k] - r_in),
                         vm_max_th_deg=float(math.degrees(geo["thc"][k])), vm_max_off_channel=float(d["vm"][k2]),
                         vm_hot_max=float(d["vm"][inn].max()), stt_mid=float(d["stt"][mid].mean()),
                         szz_mid=float(d["szz"][mid].mean()), vm_mid=float(d["vm"][mid].mean()),
                         stt_web=float(d["stt"][web].mean()), szz_web=float(d["szz"][web].mean()),
                         vm_web=float(d["vm"][web].mean()),
                         stt_jacket=float(np.average(d["stt"][jack], weights=geo["area"][jack])) if jack.any() else 0.0,
                         szz_jacket=float(np.average(d["szz"][jack], weights=geo["area"][jack])) if jack.any() else 0.0,
                         eps_z=float(d["ez"]))
    if keep_field:
        out["field"] = dict(XN=M["XN"], YN=M["YN"], conn=M["conn"], T=T, vm=st["tot"]["vm"], stt=st["tot"]["stt"],
                            szz=st["tot"]["szz"], r_in=r_in, apex=M["apex"], t=t)
    return out
