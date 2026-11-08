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
"""
import csv
import math

import numpy as np

from .channel import geometry, section
from .correlations import coolant_nu_f, curv_factor

KEYS = ("Tb", "P", "Twg", "Twc", "Taw", "pg", "q", "hg", "hc", "Re", "v", "Dh", "Fc", "eta", "rho", "reg", "ns",
        "we", "web", "sig")


class ChannelClosed(ValueError):
    pass


def solve(D, contour, gas, fluid, mdot, Tin, Pin, ks=30e-6, K_in=0.3, K_out=1.5, C_hg=1.0, k_wall=320.0,
          dx=0.25, G=None, GG=None):
    """Counter-flow march.

    D: ChannelDesign, contour: Contour, gas: IdealGas/CanteraGas, fluid: TabulatedFluid,
    mdot: total coolant flow [kg/s], Tin [K], Pin: static pressure in the inlet manifold [Pa], ks: roughness [m].
    G (geometry) and GG (gas.on_grid) can be passed in to save time in parameter studies.
    Arrays are indexed along G['x'] (injector -> exit). P is the channel static pressure [Pa].
    """
    if G is None:
        G = geometry(D, contour, dx)
    xg = G["x"]
    n = len(xg)
    dxg = float(xg[1] - xg[0])
    if GG is None:
        GG = gas.on_grid(xg)
    r, rp, ca, t, ds, kn, kg = (G[k] for k in ("r", "rp", "ca", "t", "ds", "kn", "kg"))
    N, nch = D.Neff, D.nch
    mch = mdot / nch
    Ah = 2 * math.pi * r * 1e-3 / N * np.sqrt(1 + rp ** 2) * 1e-3     # hot-wall area per pitch per mm of x [m^2]
    out = {k: np.zeros(n) for k in KEYS}
    H = fluid.H_of_T(Tin)
    P = Pin
    Q = L = dP_inl = 0.0
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
        Twg, Twc = 700.0, T + 60
        for _ in range(60):
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
            Rw = S["teff"] * 1e-3 / (k_wall * Am)
            Rc = 1 / (hc * Aw)
            Rg = 1 / (hg * Ah[j])
            qp = (Taw - T) / (Rg + Rw + Rc)
            Twg_n, Twc_n = Taw - qp * Rg, T + qp * Rc
            if abs(Twg_n - Twg) < 0.05 and abs(Twc_n - Twc) < 0.05:
                Twg, Twc = Twg_n, Twc_n
                break
            Twg = 0.6 * Twg + 0.4 * Twg_n
            Twc = 0.6 * Twc + 0.4 * Twc_n
        Ps = P - dP_inl
        pg = GG["p"][j]
        dpw = max(Ps - pg, 0.0)
        sig = dpw * (S["we"] / t[j]) ** 2 / 2 / 1e6 if D.shape == "rect" else 0.0   # hot-wall plate bending [MPa]
        vals = (T, Ps, Twg, Twc, Taw, pg, qp / Ah[j], hg, hc, Re, v, Dh * 1e3, Fc, eta, rho, 0 if reg == "lam" else 1,
                ns, S["we"], S["web"], sig)
        for k_, v_ in zip(KEYS, vals):
            out[k_][j] = v_
        if j == 0:
            break
        Q += qp * N * dxg
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
    res.update(mdot=mdot, Tin=Tin, Pin=Pin, Q_kW=Q / 1e3, dP_inl=dP_inl / 1e5, dP_out=dP_out / 1e5,
               dP_ch=(Pin - dP_inl - out["P"][0]) / 1e5,
               dP_tot=(Pin - out["P"][0] + dP_out) / 1e5,
               Pmin=float(out["P"].min()) / 1e5,
               Tout=out["Tb"][0], dh_coolant=fluid.H_of_T(out["Tb"][0]) - fluid.H_of_T(Tin), L_path=L,
               C_hg=C_hg, k_wall=k_wall)
    return res


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
    if fluid is not None and fluid.Pc:
        m["P_over_Pc"] = m["Pmin"] - fluid.Pc / 1e5
        if None not in (fluid.Tc, fluid.omega):
            sub = res["Tb"] < fluid.Tc
            if sub.any():
                m["psat_margin"] = float(np.min(res["P"][sub] / 1e5 - fluid.psat_bar(res["Tb"][sub])))
    return m


def to_csv(res, path):
    G = res["G"]
    cols = [("x_mm", res["x"]), ("r_hot_mm", G["r"]), ("t_mm", G["t"]), ("s_mm", G["s"]), ("beta_deg", G["beta"]),
            ("pitch_mm", G["pitch"]), ("web_mm", res["web"]), ("n_sub", res["ns"]), ("Tb_K", res["Tb"]),
            ("P_bar", res["P"] / 1e5), ("p_gas_bar", res["pg"] / 1e5), ("Taw_K", res["Taw"]), ("Twg_K", res["Twg"]),
            ("Twc_K", res["Twc"]), ("q_MW_m2", res["q"] / 1e6), ("hg_W_m2K", res["hg"]), ("hc_W_m2K", res["hc"]),
            ("Re", res["Re"]), ("v_m_s", res["v"]), ("Dh_mm", res["Dh"]), ("eta_fin", res["eta"]),
            ("Dean_factor", res["Fc"]), ("turbulent", res["reg"])]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([c for c, _ in cols])
        for i in range(len(res["x"])):
            w.writerow(["%.6g" % a[i] for _, a in cols])
