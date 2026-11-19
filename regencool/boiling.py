"""Boiling margins of a single-phase 1D solution (post-processing, the heat transfer itself stays single phase).

Wherever the channel static pressure is below the coolant critical pressure, the coolant-side wall (channel floor,
Twc) is compared with the saturation state at the local pressure:

  Tsat margin   Tsat(P) - Twc                                     (> 0: no vapour can form at the wall)
  ONB           Davis & Anderson (1966) wall superheat for the onset of nucleate boiling
                    dT_ONB = sqrt(8 sigma Tsat v_fg q_c / (k_l h_fg))                     uncertainty ~ +-40 %
                r_onb = (Twc - Tsat) / dT_ONB   (>= 1: nucleate boiling starts)
  CHF (best)    Hall & Mudawar (2000) subcooled flow-boiling CHF, tube inlet-conditions form
                    q_CHF = G h_fg 0.0722 We^-0.312 (rho_f/rho_g)^-0.644 [1 - 0.900 (rho_f/rho_g)^0.724 x_e]
                    We = G^2 D_h / (rho_f sigma),  x_e = -cp_l (Tsat - Tb) / h_fg  (< 0 subcooled)
                (water database, 2-200 bar; other fluids are an extrapolation, ~ +-50 %)
  CHF (bound)   Zuber pool-boiling CHF with the Ivey & Morris subcooling factor (no flow benefit -> lower bound)
                    q_Z = 0.131 rho_g^0.5 h_fg [sigma g (rho_f - rho_g)]^0.25
                          [1 + 0.1 (rho_f/rho_g)^0.75 cp_l (Tsat - Tb) / h_fg]                 ~ +-35 %
  q_c           coolant-side wall flux hc (Twc - Tb) at the channel floor.

flag: 0 single phase; 1 ONB reached; 2 CHF approach = wall at or above Tsat AND (q_c >= 0.5 q_CHF,best OR
q_c >= q_CHF,bound). A wall below Tsat generates no vapour, so a high q_c / q_CHF there is not a burnout risk.

Saturation properties come from a Saturation table: built in for isopropanol (ipa_saturation), or a CSV
(Saturation.from_csv) for any other coolant. Stations at P >= Pc, or with h_fg < 1 kJ/kg, are skipped (NaN).
"""
import csv
import math

import numpy as np

from . import fluids

G0 = 9.80665
COLS = ("P_Pa", "Tsat_K", "hfg_J_kg", "rho_l", "rho_v", "k_l", "cp_l", "sigma_N_m")


class Saturation:
    """Saturation properties tabulated against pressure (interpolated in ln P)."""

    def __init__(self, P, Tsat, hfg, rho_l, rho_v, k_l, cp_l, sigma, Pc=None, name=""):
        o = np.argsort(P)
        self.lnP = np.log(np.asarray(P, float)[o])
        self.tab = {k: np.asarray(v, float)[o] for k, v in
                    dict(Tsat=Tsat, hfg=hfg, rho_l=rho_l, rho_v=rho_v, k_l=k_l, cp_l=cp_l, sigma=sigma).items()}
        self.Pc = Pc if Pc is not None else float(np.exp(self.lnP[-1]))
        self.name = name

    def at(self, P):
        lp = math.log(P)
        if lp < self.lnP[0] or lp > self.lnP[-1]:
            return None
        return {k: float(np.interp(lp, self.lnP, v)) for k, v in self.tab.items()}

    @classmethod
    def from_csv(cls, path, Pc=None, name="csv"):
        """Columns P_bar, Tsat_K, hfg_J_kg, rho_l, rho_v, k_l, cp_l, sigma_N_m (one row per pressure)."""
        with open(path, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        col = lambda c: np.array([float(o[c]) for o in rows])  # noqa: E731
        return cls(col("P_bar") * 1e5, col("Tsat_K"), col("hfg_J_kg"), col("rho_l"), col("rho_v"), col("k_l"),
                   col("cp_l"), col("sigma_N_m"), Pc=Pc, name=name)


# ---------------------------------------------------------------------------------------------- isopropanol
IPA_TNB, IPA_DHV_NB = 355.4, 39.85e3      # normal boiling point [K], latent heat there [J/mol]
IPA_SIGMA_293 = 21.7e-3                   # surface tension at 293.15 K [N/m]


def ipa_sigma(T):
    """Surface tension, corresponding-states scaling sigma ~ (Tc - T)^(11/9) through the 293 K value."""
    return IPA_SIGMA_293 * (max(fluids.IPA_TC - T, 0.0) / (fluids.IPA_TC - 293.15)) ** (11 / 9)


def ipa_hfg(T):
    """Watson latent heat [J/kg]."""
    return IPA_DHV_NB * ((1 - T / fluids.IPA_TC) / (1 - IPA_TNB / fluids.IPA_TC)) ** 0.38 / (fluids.IPA_MW * 1e-3)


def _pr_vapour_density(T, P):
    """Peng-Robinson vapour root (largest Z) for IPA; ideal gas if only one real root."""
    om, Tc, Pc, R = fluids.IPA_OM, fluids.IPA_TC, fluids.IPA_PC, fluids.R_U
    kap = 0.379642 + 1.48503 * om - 0.164423 * om ** 2 + 0.016666 * om ** 3
    a = 0.45724 * R * R * Tc * Tc / Pc * (1 + kap * (1 - math.sqrt(T / Tc))) ** 2
    b = 0.07780 * R * Tc / Pc
    A = a * P / (R * T) ** 2
    B = b * P / (R * T)
    rr = [z.real for z in np.roots([1, -(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)])
          if abs(z.imag) < 1e-9 and z.real > B]
    Z = max(rr) if rr else 1.0
    return P * fluids.IPA_MW * 1e-3 / (Z * R * T)


def ipa_saturation(P_min=0.5e5, n=160):
    """IPA saturation table from P_min up to 0.995 Pc: Lee-Kesler Tsat, Watson h_fg, DIPPR-type liquid
    density / conductivity / cp, Peng-Robinson vapour density, scaled surface tension."""
    Tc = fluids.IPA_TC
    T = np.linspace(260.0, Tc - 0.5, 600)
    tmp = fluids.TabulatedFluid("ipa", [0, 1], [0, 1], [1, 1], [1, 1], [1, 1], Tc=Tc, Pc=fluids.IPA_PC,
                                omega=fluids.IPA_OM)
    Ps = tmp.psat_bar(T) * 1e5
    P = np.exp(np.linspace(math.log(max(P_min, Ps[0])), math.log(0.995 * fluids.IPA_PC), n))
    Ts = np.interp(np.log(P), np.log(Ps), T)
    rows = [(p, t, ipa_hfg(t), fluids._ipa_rho_liq(t), _pr_vapour_density(t, p), fluids._ipa_k_liq(t),
             fluids._ipa_cp_liq(t), ipa_sigma(t)) for p, t in zip(P, Ts)]
    c = list(zip(*rows))
    return Saturation(*c, Pc=fluids.IPA_PC, name="isopropanol")


def make(cfg, fluid, base_path=lambda p: p):
    """Saturation table for the [coolant] config: built-in IPA, or `saturation_table` CSV, else None."""
    if "saturation_table" in cfg:
        return Saturation.from_csv(base_path(cfg["saturation_table"]), Pc=fluid.Pc, name=fluid.name)
    if cfg.get("fluid", "isopropanol").lower() in ("ipa", "isopropanol", "2-propanol"):
        return ipa_saturation()
    return None


# ---------------------------------------------------------------------------------------------- check
def check(res, sat):
    """Boiling margins along a 1D solution. Returns dict of arrays (NaN where skipped) + flag (int)."""
    P, Tb, Twc, hc = res["P"], res["Tb"], res["Twc"], res["hc"]
    G = res["rho"] * res["v"]
    Dh = res["Dh"] * 1e-3
    n = len(P)
    keys = ("Tsat", "Tsat_margin", "q_c", "dT_onb", "r_onb", "x_e", "q_chf", "q_chf_bound", "r_chf", "r_chf_bound")
    out = {k: np.full(n, np.nan) for k in keys}
    flag = np.zeros(n, dtype=int)
    for j in range(n):
        if P[j] >= sat.Pc:
            continue
        s = sat.at(P[j])
        if s is None or s["hfg"] < 1e3:
            continue
        Ts, hfg, rf, rg, kl, cpl, sg = (s[k] for k in ("Tsat", "hfg", "rho_l", "rho_v", "k_l", "cp_l", "sigma"))
        qc = max(hc[j] * (Twc[j] - Tb[j]), 0.0)
        dTsub = Ts - Tb[j]
        xe = -cpl * dTsub / hfg
        vfg = 1 / rg - 1 / rf
        dT_onb = math.sqrt(8 * sg * Ts * vfg * max(qc, 1.0) / (kl * hfg)) if sg > 0 else 0.0
        dr = rf / rg
        We = G[j] ** 2 * Dh[j] / (rf * sg) if sg > 0 else np.inf
        q_hm = G[j] * hfg * 0.0722 * We ** -0.312 * dr ** -0.644 * (1 - 0.900 * dr ** 0.724 * min(xe, 0.0))
        q_z = 0.131 * rg ** 0.5 * hfg * (sg * G0 * (rf - rg)) ** 0.25 * (1 + 0.1 * dr ** 0.75 * cpl * max(dTsub, 0)
                                                                        / hfg)
        vals = dict(Tsat=Ts, Tsat_margin=Ts - Twc[j], q_c=qc, dT_onb=dT_onb,
                    r_onb=(Twc[j] - Ts) / dT_onb if dT_onb > 0 else np.inf, x_e=xe, q_chf=q_hm, q_chf_bound=q_z,
                    r_chf=qc / q_hm if q_hm > 0 else np.inf, r_chf_bound=qc / q_z if q_z > 0 else np.inf)
        for k, v in vals.items():
            out[k][j] = v
        if vals["r_onb"] >= 1.0:
            flag[j] = 1
        if Twc[j] >= Ts and (vals["r_chf"] >= 0.5 or vals["r_chf_bound"] >= 1.0):
            flag[j] = 2
    out["flag"] = flag
    return out


def summary(b, x):
    """Scalar summary of check(): worst values and where."""
    sub = np.isfinite(b["Tsat_margin"])
    if not sub.any():
        return dict(subcritical_stations=0, flag_max=0)
    i = int(np.nanargmin(b["Tsat_margin"]))
    m = dict(subcritical_stations=int(sub.sum()), Tsat_margin_min=float(b["Tsat_margin"][i]), x_Tsat_margin=float(x[i]),
             r_onb_max=float(np.nanmax(b["r_onb"])), r_chf_max=float(np.nanmax(b["r_chf"])),
             r_chf_bound_max=float(np.nanmax(b["r_chf_bound"])), flag_max=int(b["flag"].max()))
    for f in (1, 2):
        hit = np.where(b["flag"] >= f)[0]
        if len(hit):
            m["x_flag%d" % f] = [float(x[hit].min()), float(x[hit].max())]
    return m


def csv_columns(b):
    return [("Tsat_K", b["Tsat"]), ("Tsat_margin_K", b["Tsat_margin"]), ("q_c_MW_m2", b["q_c"] / 1e6),
            ("r_onb", b["r_onb"]), ("q_chf_MW_m2", b["q_chf"] / 1e6), ("r_chf", b["r_chf"]),
            ("r_chf_bound", b["r_chf_bound"]), ("boil_flag", b["flag"].astype(float))]
