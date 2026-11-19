"""Order-of-magnitude estimate of hot-gas radiation to the wall (post-processing of a 1D solution).

Gas emissivity: Leckner (1972) total-emissivity correlations for H2O and CO2 with pressure correction and H2O/CO2
band-overlap correction, as given in Modest, Radiative Heat Transfer. Leckner is fitted for 400-2500 K;
combustion gases are hotter, so the result is reported with T clipped at T_clip (default 2500 K) and, if
requested, extrapolated - the spread between the two is part of the answer.

Along the contour: gas at the local static temperature and pressure, mean beam length L_m = 0.9 D (cylinder),
gas absorptivity by Hottel's rule a_g = eps_g(T_w, p_a L T_w/T_g) (T_g/T_w)^0.5, grey wall with eps_w:

    q_rad = (1 + eps_w)/2 * sigma (eps_g T_g^4 - a_g T_w^4)

The convective solution is not re-solved. The wall heats up by dT = q_rad R_c / (1 + h_g R_c) with
R_c = (Twg - Tb) / q (wall + coolant resistance per area from the 1D result), i.e. the extra flux is split between
the gas side (smaller Taw - Twg) and the coolant side. Use it to see whether radiation matters (typically a few
percent of the convective flux in small engines) - not as a design value.
"""
import math

import numpy as np

SIG = 5.670374e-8

_CH2O = np.array([[-2.2118, -1.1987, 0.035596],
                  [0.85667, 0.93048, -0.14391],
                  [-0.10838, -0.17156, 0.045915]])         # rows: power of log10(p_a L), cols: power of T/1000
_CCO2 = np.array([[-3.9893, 2.7669, -2.1081, 0.39163],
                  [1.2710, -1.1090, 1.0195, -0.21897],
                  [-0.23678, 0.19731, -0.19544, 0.044644]])


def _eps0(C, paL, t):
    lg = math.log10(paL)
    a = [np.polyval(C[i][::-1], t) for i in range(C.shape[0])]
    return math.exp(a[0] + sum(a[i] * lg ** i for i in range(1, len(a))))


def eps_h2o(pa, p, L, T):
    """H2O emissivity. pa, p [bar], L [cm], T [K]."""
    t = T / 1000.0
    paL = max(pa * L, 1e-6)
    PE = p + 2.56 * pa / math.sqrt(t)
    a = 1.888 - 2.053 * math.log10(t) if t > 0.75 else 2.144
    b, c, paLm = 1.10 / t ** 1.4, 0.5, 13.2 * t * t
    corr = 1 - (a - 1) * (1 - PE) / (a + b - 1 + PE) * math.exp(-c * math.log10(paLm / paL) ** 2)
    return _eps0(_CH2O, paL, t) * corr


def eps_co2(pa, p, L, T):
    """CO2 emissivity. pa, p [bar], L [cm], T [K]."""
    t = T / 1000.0
    paL = max(pa * L, 1e-6)
    PE = p + 0.28 * pa
    a, b, c = 1 + 0.1 / t ** 1.45, 0.23, 1.47
    paLm = 0.054 / t ** 2 if t < 0.7 else 0.225 * t * t
    corr = 1 - (a - 1) * (1 - PE) / (a + b - 1 + PE) * math.exp(-c * math.log10(paLm / paL) ** 2)
    return _eps0(_CCO2, paL, t) * corr


def eps_mix(x_h2o, x_co2, p, L, T):
    """Mixture emissivity with the overlap correction. Returns (eps, eps_H2O, eps_CO2, overlap)."""
    pw, pc = x_h2o * p, x_co2 * p
    ew = eps_h2o(pw, p, L, T) if pw > 0 else 0.0
    ec = eps_co2(pc, p, L, T) if pc > 0 else 0.0
    if pw + pc <= 0:
        return 0.0, 0.0, 0.0, 0.0
    z = pw / (pw + pc)
    lg = math.log10((pw + pc) * L)
    de = (z / (10.7 + 101 * z) - 0.0089 * z ** 10.4) * max(lg, 0.0) ** 2.76
    return ew + ec - de, ew, ec, de


def estimate(res, GG, x_h2o, x_co2, eps_wall=0.8, T_clip=2500.0, beam=0.9):
    """Radiation along a 1D solution. GG: gas.on_grid (needs T, p). Returns dict of arrays + totals."""
    x, r = res["x"], res["G"]["r"]
    Tg, p = np.asarray(GG["T"], float), np.asarray(GG["p"], float) / 1e5
    Twg, Tb, q, hg = res["Twg"], res["Tb"], res["q"], res["hg"]
    Te = np.minimum(Tg, T_clip) if T_clip else Tg
    Lm = beam * 2 * r / 10.0                                            # cm
    eg = np.array([eps_mix(x_h2o, x_co2, pi, Li, Ti)[0] for pi, Li, Ti in zip(p, Lm, Te)])
    ag = np.array([eps_mix(x_h2o, x_co2, pi, Li * tw / tg, min(tw, T_clip or tw))[0] * math.sqrt(tg / tw)
                   for pi, Li, tw, tg in zip(p, Lm, Twg, Tg)])
    ag = np.clip(ag, 0.0, 1.0)
    qr = (1 + eps_wall) / 2 * SIG * (eg * Tg ** 4 - ag * Twg ** 4)
    Rc = (Twg - Tb) / np.maximum(q, 1.0)
    dT = qr * Rc / (1 + hg * Rc)
    ds = np.sqrt(1 + np.gradient(r, x) ** 2) * np.gradient(x) * 1e-3
    Q = float(np.sum(qr * 2 * math.pi * r * 1e-3 * ds))
    i = int(np.argmax(Twg + dT))
    return dict(eps_g=eg, a_g=ag, q_rad=qr, dTwg=dT, Q_rad_kW=Q / 1e3, ratio_max=float(np.max(qr / np.maximum(q, 1))),
                q_rad_max=float(qr.max()), dTwg_max=float(dT.max()), x_dTwg_max=float(x[int(np.argmax(dT))]),
                Twg_max_with_rad=float((Twg + dT)[i] - 273.15), x_Twg_max_with_rad=float(x[i]))
