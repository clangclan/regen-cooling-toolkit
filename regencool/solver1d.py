"""1D counter-flow regenerative-cooling march.

The coolant enters at the nozzle-exit end (largest x) and leaves at the injector end (x0). At each station the
heat flux per channel pitch follows from three resistances in series

    q' = (Taw - Tb) / (R_gas + R_wall + R_cool)
    R_gas  = 1 / (hg A_hot)                     hg = C_hg * Bartz h0 * sigma(Twg)
    R_wall = t_eff / (k_wall A_mid)
    R_cool = 1 / (hc (P_direct + eta_fin P_fin) ds)

solved by fixed-point iteration on Twg / Twc (relaxation 0.6/0.4). The coolant is marched in enthalpy, so the
liquid -> pseudo-critical -> supercritical transition is continuous. Pressure: Darcy friction with a Dean factor on
curved paths + acceleration (momentum) term; inlet/outlet minor losses K_in, K_out on the channel dynamic pressure.

Optional axial wall conduction (axial=True): after each march the wall temperatures of all stations are solved
together as one linear conduction network (the march's frozen gas / wall / coolant conductances + axial conduction
in the hot-wall slab and in the solid between the channels). The march is then repeated with the axial terms
implicit in each station's own temperatures and the neighbours frozen at the relaxed global solution (Aitken
dynamic relaxation) until march and global system agree. See docs/2_regen_1d.md, section 2.8.
"""
import csv
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .channel import geometry, section
from .correlations import coolant_nu_f, curv_factor

KEYS = ("Tb", "P", "Twg", "Twc", "Taw", "pg", "q", "hg", "hc", "Re", "v", "Dh", "Fc", "eta", "rho", "reg", "ns",
        "we", "web", "sig", "qc", "Ug", "Uw", "Uc", "Asol")
AX_COEFS = ("g0", "ggg", "ggc", "c0", "cgg", "ccc")
M_BIL = ((1 / 3, 1 / 6), (1 / 6, 1 / 3))   # through-thickness weights of a bilinear wall element


class ChannelClosed(ValueError):
    pass


def solve(D, contour, gas, fluid, mdot, Tin, Pin, ks=30e-6, K_in=0.3, K_out=1.5, C_hg=1.0, k_wall=320.0,
          dx=0.25, G=None, GG=None, axial=False, ax_tol=0.3, ax_maxit=30):
    """Counter-flow march.

    D: ChannelDesign, contour: Contour, gas: IdealGas/CanteraGas, fluid: TabulatedFluid,
    mdot: total coolant flow [kg/s], Tin [K], Pin: static pressure in the inlet manifold [Pa], ks: roughness [m].
    G (geometry) and GG (gas.on_grid) can be passed in to save time in parameter studies.
    axial: couple axial wall conduction (iterate until max |Twg(march) - Twg(global)| < ax_tol [K]).
    Arrays are indexed along G['x'] (injector -> exit). P is the channel static pressure [Pa].
    """
    if G is None:
        G = geometry(D, contour, dx)
    if GG is None:
        GG = gas.on_grid(G["x"])
    args = (D, G, GG, gas, fluid, mdot, Tin, Pin, ks, K_in, K_out, C_hg, k_wall)
    res = _march(*args)
    if not axial:
        return res
    n = len(G["x"])
    dxg = float(G["x"][1] - G["x"][0])
    first = res
    hist = []
    xk = r_prev = None
    w = 1.0
    conv = False
    for it in range(ax_maxit):
        els = axial_elements(D, G, res, k_wall)
        Tg_s, Tc_s = axial_system(G, res, els)
        dsys = float(np.max(np.abs(res["Twg"] - Tg_s)))
        hist.append(dict(it=it, dsys=dsys, relax=w, Q_kW=res["Q_kW"], Qg_kW=res["Qg_kW"]))
        if dsys < ax_tol:
            conv = True
            break
        xs = np.concatenate([Tg_s, Tc_s])
        if xk is None:
            xk = xs
        else:
            rk = xs - xk
            if r_prev is not None:                       # Aitken dynamic relaxation (Irons & Tuck)
                dr = rk - r_prev
                den = float(dr @ dr)
                if den > 0:
                    w = min(max(-w * float(r_prev @ dr) / den, 0.1), 1.5)
            xk = xk + w * rk
            r_prev = rk
        XC = axial_coefs(n, els, xk[:n], xk[n:], res["Tb"], dxg)
        res = _march(*args, XC=XC)
    res.update(axial=True, axial_conv=conv, axial_hist=hist, Twg_noax=first["Twg"], Twc_noax=first["Twc"])
    return res


def _march(D, G, GG, gas, fluid, mdot, Tin, Pin, ks, K_in, K_out, C_hg, k_wall, XC=None):
    """One counter-flow march. XC: implicit axial-conduction coefficients per station (axial_coefs) or None."""
    xg = G["x"]
    n = len(xg)
    dxg = float(xg[1] - xg[0])
    r, rp, ca, t, ds, kn, kg = (G[k] for k in ("r", "rp", "ca", "t", "ds", "kn", "kg"))
    N, nch = D.Neff, D.nch
    mch = mdot / nch
    Ah = 2 * math.pi * r * 1e-3 / N * np.sqrt(1 + rp ** 2) * 1e-3     # hot-wall area per pitch per mm of x [m^2]
    out = {k: np.zeros(n) for k in KEYS}
    H = fluid.H_of_T(Tin)
    P = Pin
    Q = Qg = L = dP_inl = 0.0
    zero = dict.fromkeys(AX_COEFS, 0.0)
    for j in range(n - 1, -1, -1):
        T = fluid.T_of_H(H)
        S = section(D, G, j)
        ns = S["ns"]
        mce = mch / ns
        if S["we"] <= 0.05 or S["hh"] <= 0.05:
            raise ChannelClosed("channel closed at x = %.2f mm" % xg[j])
        Dh = S["Dh"] * 1e-3
        A = S["A"] * 1e-6
        rho, mu, cp, kf = (fluid.at(T, k) for k in ("rho", "mu", "cp", "k"))
        v = mce / (rho * A)
        Re = rho * v * Dh / mu
        Pr = mu * cp / kf
        if j == n - 1:
            dP_inl = K_in * 0.5 * rho * v * v
        h0, Taw, g1 = GG["h0"][j], GG["Taw"][j], GG["g1"][j]
        X = zero if XC is None else {k_: float(XC[k_][j]) for k_ in AX_COEFS}
        Twg, Twc = 700.0, T + 60
        for _ in range(80):
            hg = C_hg * h0 * gas.sigma(Twg, g1)
            f, Nu, reg = coolant_nu_f(fluid, Re, Pr, Dh, S["a"], ks, L * 1e-3 + 1e-3, T, max(Twc, T + 0.1), mu,
                                      fluid.at(max(Twc, T), "mu"))
            if D.shape != "rect" and reg == "lam":       # circular-duct laminar constants
                f *= 64 / 56.9
                Nu *= 4.36 / 3.61
            Fc = curv_factor(Re, Dh, kn[j] * 1e3) if reg == "turb" else 1.0
            hc = Nu * kf / Dh * Fc
            m = math.sqrt(2 * hc / (k_wall * S["tf"] * 1e-3))
            mh = m * S["Hf"] * 1e-3
            eta = math.tanh(mh) / mh
            Aw = ns * (S["Pd"] + eta * S["Pf"]) * 1e-3 * ds[j] * 1e-3
            Am = 2 * math.pi * (r[j] + 0.5 * t[j] * ca[j]) * 1e-3 / N * math.sqrt(1 + rp[j] ** 2) * 1e-3
            a_ = hg * Ah[j]                               # gas, wall, coolant conductance [W/K per mm of x]
            b_ = k_wall * Am / (S["teff"] * 1e-3)
            c_ = hc * Aw
            # 2x2 station network; with XC the axial heat from the frozen neighbours enters implicitly
            A11, A12 = a_ + b_ + X["ggg"], -b_ + X["ggc"]
            A21, A22 = -b_ + X["cgg"], b_ + c_ + X["ccc"]
            r1, r2 = a_ * Taw + X["g0"], c_ * T + X["c0"]
            det = A11 * A22 - A12 * A21
            Twg_n = (r1 * A22 - A12 * r2) / det
            Twc_n = (A11 * r2 - A21 * r1) / det
            if abs(Twg_n - Twg) < 0.05 and abs(Twc_n - Twc) < 0.05:
                Twg, Twc = Twg_n, Twc_n
                break
            Twg = 0.6 * Twg + 0.4 * Twg_n
            Twc = 0.6 * Twc + 0.4 * Twc_n
        qp = c_ * (Twc - T)                               # heat into the coolant [W per mm of x per pitch]
        qg = a_ * (Taw - Twg)                             # heat from the gas (differs from qp by the axial heat)
        Ps = P - dP_inl
        pg = GG["p"][j]
        dpw = max(Ps - pg, 0.0)
        sig = dpw * (S["we"] / t[j]) ** 2 / 2 / 1e6 if D.shape == "rect" else 0.0   # hot-wall plate bending [MPa]
        Asol = max(G["pitch"][j] * S["hh"] - ns * S["A"], 0.0)   # solid between the channels [mm^2 per pitch]
        vals = (T, Ps, Twg, Twc, Taw, pg, qg / Ah[j], hg, hc, Re, v, Dh * 1e3, Fc, eta, rho, 0 if reg == "lam" else 1,
                ns, S["we"], S["web"], sig, qp / Ah[j], a_, b_, c_, Asol)
        for k_, v_ in zip(KEYS, vals):
            out[k_][j] = v_
        if j == 0:
            break
        Q += qp * N * dxg
        Qg += qg * N * dxg
        H += qp * N / nch * dxg / mch
        kc = math.sqrt(kn[j] ** 2 + kg[j] ** 2) * 1e3
        fcf = 1.0
        if reg == "turb" and kc > 0:
            Dn = Re * (Dh * kc / 2) ** 2
            if Dn > 1:
                fcf = Dn ** 0.05
        dPf = f * fcf * (ds[j] * dxg * 1e-3) / Dh * 0.5 * rho * v * v
        rho_n = fluid.at(fluid.T_of_H(H), "rho")
        dPa = (mce / A) ** 2 * (1 / rho_n - 1 / rho)
        P -= dPf + dPa
        L += ds[j] * dxg
    dP_out = K_out * 0.5 * out["rho"][0] * out["v"][0] ** 2
    res = dict(G=G, x=xg, **out)
    res.update(mdot=mdot, Tin=Tin, Pin=Pin, Q_kW=Q / 1e3, Qg_kW=Qg / 1e3, axial=False, dP_inl=dP_inl / 1e5,
               dP_out=dP_out / 1e5,
               dP_ch=(Pin - dP_inl - out["P"][0]) / 1e5,
               dP_tot=(Pin - out["P"][0] + dP_out) / 1e5,
               Pmin=float(out["P"].min()) / 1e5,
               Tout=out["Tb"][0], dh_coolant=fluid.H_of_T(out["Tb"][0]) - fluid.H_of_T(Tin), L_path=L,
               C_hg=C_hg, k_wall=k_wall)
    return res


# ------------------------------------------------------------------------------------- axial wall conduction
def _mid_length(G):
    """Meridional length of the wall mid-surface per unit x (contour offset by t/2 along its normal)."""
    x, r, rp, t = G["x"], G["r"], G["rp"], G["t"]
    ph = np.arctan(rp)
    xm = x - 0.5 * t * np.sin(ph)
    rm = r + 0.5 * t * np.cos(ph)
    return np.hypot(np.gradient(xm, x), np.gradient(rm, x))


def axial_elements(D, G, res, k_wall):
    """Axial conductances [W/K] between neighbouring stations j, j+1, per channel pitch.

    Ks: hot-wall slab (thickness t, width and length on the wall mid-surface).
    Kr: solid between the channels (ribs, lands, webs); it conducts along the channel path, so its meridional
        conductance carries a factor cos(beta) on helical channels.
    eta: mean fin efficiency of the two stations (mean temperature excess of the solid over the bulk, relative
        to the channel floor)."""
    x = G["x"]
    dx = float(x[1] - x[0])
    r, t, ca, cosb, ds = G["r"], G["t"], G["ca"], G["cosb"], G["ds"]
    wmid = 2 * math.pi * (r + 0.5 * t * ca) / D.Neff                  # mm
    smid = _mid_length(G)
    As = wmid * t                                                      # mm^2
    Ar = res["Asol"]
    els = []
    for j in range(len(x) - 1):
        k = j + 1
        Le = 0.5 * (smid[j] + smid[k]) * dx
        Ll = 0.5 * (ds[j] * cosb[j] + ds[k] * cosb[k]) * dx
        cb = 0.5 * (cosb[j] + cosb[k])
        Ks = k_wall * 0.5 * (As[j] + As[k]) * 1e-3 / Le                  # k A / L with A [mm^2], L [mm]
        Kr = k_wall * 0.5 * (Ar[j] + Ar[k]) * cb * 1e-3 / Ll
        eta = 0.5 * (res["eta"][j] + res["eta"][k])
        els.append((j, k, Ks, Kr, eta))
    return els


def axial_system(G, res, els):
    """Global wall network with the march's conductances frozen. Unknowns: (Twg_j, Twc_j) of all stations.

    Station balance per pitch [W]: gas a (Taw - Twg), wall b (Twc - Twg), coolant c (Tb - Twc), plus the axial heat
    from the neighbours. Slab: bilinear element (weights 1/3, 1/6 through the thickness). Solid between the
    channels: mean temperature Tb + eta (Twc - Tb) (fin profile); its axial heat is booked on the coolant-side node.
    Ends adiabatic. Returns (Twg, Twc)."""
    n = len(G["x"])
    dx = float(G["x"][1] - G["x"][0])
    a, b, c = res["Ug"] * dx, res["Uw"] * dx, res["Uc"] * dx
    Taw, Tb = res["Taw"], res["Tb"]
    gi = 2 * np.arange(n)
    ci = gi + 1
    rows, cols, vals = [], [], []

    def add(i, j_, v):
        rows.append(i)
        cols.append(j_)
        vals.append(v)

    rhs = np.zeros(2 * n)
    for p in range(n):
        add(gi[p], gi[p], -(a[p] + b[p]))
        add(gi[p], ci[p], b[p])
        add(ci[p], gi[p], b[p])
        add(ci[p], ci[p], -(b[p] + c[p]))
        rhs[gi[p]] = -a[p] * Taw[p]
        rhs[ci[p]] = -c[p] * Tb[p]
    for j, k, Ks, Kr, eta in els:
        for p, q in ((j, k), (k, j)):
            for ri, rr in ((0, gi), (1, ci)):
                for cj, cc in ((0, gi), (1, ci)):
                    add(rr[p], cc[q], Ks * M_BIL[ri][cj])
                    add(rr[p], cc[p], -Ks * M_BIL[ri][cj])
            add(ci[p], ci[q], Kr * eta)
            add(ci[p], ci[p], -Kr * eta)
            rhs[ci[p]] -= Kr * (1 - eta) * (Tb[q] - Tb[p])
    A = sp.csc_matrix((vals, (rows, cols)), shape=(2 * n, 2 * n))
    T = spla.spsolve(A, rhs)
    return T[gi], T[ci]


def axial_coefs(n, els, Tg, Tc, Tb, dx):
    """Implicit axial heat into station p with the neighbours q frozen at (Tg, Tc) [W per mm of x per pitch]:
    hot-side node  g0 - ggg Twg_p - ggc Twc_p,  coolant-side node  c0 - cgg Twg_p - ccc Twc_p."""
    C = {k: np.zeros(n) for k in AX_COEFS}
    m0, m1 = M_BIL[0]
    for j, k, Ks, Kr, eta in els:
        for p, q in ((j, k), (k, j)):
            C["g0"][p] += Ks * (m0 * Tg[q] + m1 * Tc[q])
            C["ggg"][p] += Ks * m0
            C["ggc"][p] += Ks * m1
            C["c0"][p] += Ks * (m1 * Tg[q] + m0 * Tc[q]) + Kr * eta * Tc[q] + Kr * (1 - eta) * (Tb[q] - Tb[p])
            C["cgg"][p] += Ks * m1
            C["ccc"][p] += Ks * m0 + Kr * eta
    return {k: v / dx for k, v in C.items()}


def metrics(res, contour, fluid=None):
    """Scalar summary (temperatures in deg C, pressures in bar)."""
    x = res["x"]
    it = int(np.argmin(np.abs(x - contour.x_t)))
    ig = int(np.argmax(res["Twg"]))
    m = dict(dP=float(res["dP_tot"]), dP_ch=float(res["dP_ch"]), Pmin=float(res["Pmin"]),
             Twg=float(res["Twg"][ig] - 273.15), x_Twg=float(x[ig]), Twg_throat=float(res["Twg"][it] - 273.15),
             Twc=float(res["Twc"].max() - 273.15), Tout=float(res["Tout"] - 273.15), Q_kW=float(res["Q_kW"]),
             dh_kJkg=float(res["dh_coolant"] / 1e3), q_max=float(res["q"].max() / 1e6),
             v_throat=float(res["v"][it]), v_max=float(res["v"].max()), Re_min=float(res["Re"].min()),
             sig_plate=float(res["sig"].max()), web_min=float(res["web"].min()), w_max=float(res["we"].max()),
             s_min=float(res["G"]["s"].min()), beta_max=float(np.abs(res["G"]["beta"]).max()),
             laminar_frac=float(1 - res["reg"].mean()))
    if res.get("axial"):
        m["Twg_noax"] = float(res["Twg_noax"].max() - 273.15)
        m["axial_iterations"] = len(res["axial_hist"])
        m["axial_converged"] = bool(res["axial_conv"])
        m["Qg_kW"] = float(res["Qg_kW"])
    if fluid is not None and fluid.Pc:
        m["P_over_Pc"] = m["Pmin"] - fluid.Pc / 1e5
        if None not in (fluid.Tc, fluid.omega):
            sub = res["Tb"] < fluid.Tc
            if sub.any():
                m["psat_margin"] = float(np.min(res["P"][sub] / 1e5 - fluid.psat_bar(res["Tb"][sub])))
    return m


def to_csv(res, path, extra=None):
    """extra: optional list of (column name, array) appended to the standard columns (e.g. boiling checks)."""
    G = res["G"]
    cols = [("x_mm", res["x"]), ("r_hot_mm", G["r"]), ("t_mm", G["t"]), ("s_mm", G["s"]), ("beta_deg", G["beta"]),
            ("pitch_mm", G["pitch"]), ("web_mm", res["web"]), ("n_sub", res["ns"]), ("Tb_K", res["Tb"]),
            ("P_bar", res["P"] / 1e5), ("p_gas_bar", res["pg"] / 1e5), ("Taw_K", res["Taw"]), ("Twg_K", res["Twg"]),
            ("Twc_K", res["Twc"]), ("q_MW_m2", res["q"] / 1e6), ("q_cool_MW_m2", res["qc"] / 1e6),
            ("hg_W_m2K", res["hg"]), ("hc_W_m2K", res["hc"]),
            ("Re", res["Re"]), ("v_m_s", res["v"]), ("Dh_mm", res["Dh"]), ("eta_fin", res["eta"]),
            ("Dean_factor", res["Fc"]), ("turbulent", res["reg"])]
    if res.get("axial"):
        cols.append(("Twg_no_axial_K", res["Twg_noax"]))
    cols += list(extra or [])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([c for c, _ in cols])
        for i in range(len(res["x"])):
            w.writerow(["%.6g" % a[i] for _, a in cols])
