"""Command line:  python -m regencool <command> case.toml [-o OUT] [--C-hg X]

  gas        hot-gas side along the contour (+ Cantera chamber state), nozzle separation check
  run        1D regenerative-cooling solution (+ 2D thermo-elastic stations from [stress])
  optimize   channel-size distribution s(x) for minimum hot-wall temperature under dP / pressure limits
  tolerance  uniform and single-channel manufacturing offsets
  manifold   ring-header flow distribution + per-channel tolerance Monte Carlo
  radiation  hot-gas radiation estimate on top of the 1D solution
  injector   orifice sizing from [injector] (SPI, real-fluid IPA, N2O Dyer), optionally fed by the regen outlet
"""
import argparse
import csv
import json
import os
import time

import numpy as np

from . import boiling, plots, solver1d, station2d
from .case import Case, jsonable


def _fmt_metrics(m):
    s = ("  max Twg %.0f C at x = %.1f mm (throat %.0f C), max Twc %.0f C\n"
         "  coolant out %.0f C (dh %.0f kJ/kg), Q %.1f kW, q_max %.1f MW/m2\n"
         "  dP total %.2f bar (channels %.2f), channel exit P %.2f bar, v_max %.1f m/s, Re_min %.0f"
         % (m["Twg"], m["x_Twg"], m["Twg_throat"], m["Twc"], m["Tout"], m["dh_kJkg"], m["Q_kW"], m["q_max"],
            m["dP"], m["dP_ch"], m["Pmin"], m["v_max"], m["Re_min"]))
    if "P_over_Pc" in m:
        s += "\n  exit P - coolant Pc = %+.2f bar" % m["P_over_Pc"]
    if "psat_margin" in m:
        s += ", min (P - psat) = %.2f bar" % m["psat_margin"]
    if "Twg_noax" in m:
        s += "\n  axial wall conduction: %d iterations (%s), max Twg without it %.0f C" % (
            m["axial_iterations"], "converged" if m["axial_converged"] else "NOT converged", m["Twg_noax"])
    return s


def _fmt_boiling(bs):
    s = ("  boiling check (%d stations below Pc): min Tsat - Twc %.1f K at x = %.1f mm, max (Twc - Tsat)/dT_ONB %.2f, "
         "max q_c/q_CHF %.3f (bound %.3f) -> flag %d"
         % (bs["subcritical_stations"], bs["Tsat_margin_min"], bs["x_Tsat_margin"], bs["r_onb_max"], bs["r_chf_max"],
            bs["r_chf_bound_max"], bs["flag_max"]))
    for f, lab in ((1, "ONB"), (2, "CHF approach")):
        if "x_flag%d" % f in bs:
            s += "\n    %s: x = %.1f .. %.1f mm" % ((lab,) + tuple(bs["x_flag%d" % f]))
    return s


def _write_json(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(jsonable(obj), f, indent=1)


def cmd_gas(case, a):
    print(case.describe())
    if case.regen_coupled:
        case.couple()
    g = case.gas
    if hasattr(g, "table"):
        ch = g.table["chamber"]
        print("  chamber: T0 %.0f K, MW %.2f, gamma_f %.4f, gamma_s %.4f, cp_f %.0f J/kg/K, Pr_f %.3f, c* %.0f m/s"
              % (ch["T"], ch["MW"], ch["gam_f"], ch["gam_s"], ch["cp_f"], ch["Pr_f"], g.cstar_ideal))
        print("  X: " + ", ".join("%s %.3f" % kv for kv in ch["X"].items()))
    nz = case.nozzle_flow_check()
    print("  exit: p %.3f bar, M %.2f, pe/pa %.3f, Schmucker p_sep %.3f bar -> %s"
          % (nz["pe_bar"], nz["Me"], nz["pe_over_pa"], nz["p_sep_schmucker_bar"],
             "SEPARATED at sea level" if nz["separated_schmucker"] else "attached"))
    plots.gas(case, os.path.join(case.out_dir, "gas.png"))
    _write_json(os.path.join(case.out_dir, "gas_summary.json"), dict(describe=g.describe(), nozzle=nz,
                                                                      chamber=getattr(g, "table", {}).get("chamber")))


def cmd_run(case, a):
    print(case.describe())
    t0 = time.time()
    res = case.couple()
    m = case.metrics(res)
    print("1D solution (C_hg %.2f, %.1f s):" % (case.C_hg, time.time() - t0))
    print(_fmt_metrics(m))
    b = case.boiling(res)
    extra = None
    out_b = None
    if b is not None:
        out_b = boiling.summary(b, res["x"])
        print(_fmt_boiling(out_b))
        extra = boiling.csv_columns(b)
        plots.boiling(case, res, b, os.path.join(case.out_dir, "boiling.png"))
    solver1d.to_csv(res, os.path.join(case.out_dir, "solution_1d.csv"), extra)
    plots.geometry(case, res, os.path.join(case.out_dir, "geometry.png"))
    plots.solution(case, res, os.path.join(case.out_dir, "solution_1d.png"),
                   case.cfg.get("limits", {}).get("Twg_max_C"))
    out = dict(metrics=m, nozzle=case.nozzle_flow_check(), dh_fuel_kJkg=case.dh_fuel / 1e3, T0=case.gas.T0,
               boiling=out_b, axial_hist=res.get("axial_hist"))
    if not a.no_stress:
        st_cfg = case.cfg.get("stress", {})
        D = case.design()
        out["stations"] = []
        for x in case.stations(res):
            st = station2d.station(x, D, res, case.material, jacket=st_cfg.get("jacket_mm", 2.0),
                                   h=st_cfg.get("h_mm", 0.015), keep_field=True)
            print("2D station x = %.1f mm: Twg %.0f C (1D %.0f), vM max %.0f MPa (thermal %.0f, pressure %.0f), "
                  "hot-wall mid hoop %.0f MPa, axial %.0f MPa"
                  % (st["x"], st["Twg_max"], st["Twg_1d"], st["tot"]["vm_max"], st["th"]["vm_max"],
                     st["pr"]["vm_max"], st["tot"]["stt_mid"], st["tot"]["szz_mid"]))
            plots.station_field(st, os.path.join(case.out_dir, "station_x%.1f.png" % st["x"]))
            out["stations"].append(st)
        if out["stations"]:
            plots.stress_summary(out["stations"], os.path.join(case.out_dir, "stress_summary.png"))
    _write_json(os.path.join(case.out_dir, "summary.json"), out)
    print("-> %s" % case.out_dir)


def cmd_optimize(case, a):
    from .optimize import Optimizer, toml_snippet
    print(case.describe())
    base = case.couple()
    oc = case.cfg.get("optimize", {})
    opt = Optimizer(case, oc.get("knots_x_mm"), oc.get("dP_max_bar"), oc.get("P_min_bar"), oc.get("robust_d_mm", 0.0),
                    oc.get("s_min_mm", 0.3), C_hg=a.C_hg)
    print("optimizing s(x): knots %s, dP_max %.2f bar, P_min %.2f bar, robust d %.3f mm"
          % (np.round(opt.kx, 2).tolist(), opt.dP_max, opt.P_min, opt.d))
    r = opt.run(oc.get("maxfev", 200))
    print("optimized (%d evaluations):" % r["nfev"])
    print(_fmt_metrics(r["nominal"]))
    snip = toml_snippet(r["knots"])
    print("paste into [channels]:\n" + snip)
    with open(os.path.join(case.out_dir, "optimized_s.toml"), "w") as f:
        f.write("# optimized channel size knots (x_mm, s_mm)\n" + snip)
    after = case.solve(case.design(s_fun=opt.s_fun(np.log(np.array(r["knots"])[:, 1]))))
    plots.optimization(case, base, after, os.path.join(case.out_dir, "optimization.png"))
    _write_json(os.path.join(case.out_dir, "optimization.json"), r)


def cmd_tolerance(case, a):
    from .tolerance import study
    print(case.describe())
    case.couple()
    tc = case.cfg.get("tolerance", {})
    rows = study(case, tc.get("offsets_mm", [-0.2, -0.1, -0.05, 0.05, 0.1, 0.2]),
                 tc.get("one_channel_mm", [0.05, 0.1, 0.2]), C_hg=a.C_hg)
    keys = ["case", "d", "flow", "Twg", "x_Twg", "Twg_throat", "Tout", "dP", "Pmin", "v_max", "s_min", "web_min"]
    with open(os.path.join(case.out_dir, "tolerance.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] if isinstance(r[k], str) else "%.4g" % r[k] for k in keys])
    plots.tolerance(rows, os.path.join(case.out_dir, "tolerance.png"))
    _write_json(os.path.join(case.out_dir, "tolerance.json"), rows)


def cmd_manifold(case, a):
    from . import manifold as MF
    print(case.describe())
    case.couple()
    mc = case.cfg["manifold"]
    D = case.design()
    ang = MF.channel_angles(D.nch, D.NL, D.land, mc.get("land0_deg", 0.0))
    rho, mu = case.fluid.at(case.Tin, "rho"), case.fluid.at(case.Tin, "mu")
    off, prof = MF.header_offsets(ang, case.mdot, rho, mu, mc, mc.get("k_r", 0.8), mc.get("C_b", 0.5), case.ks)
    print("header: static offset %.3f .. %.3f bar over %d channels" % (off.min(), off.max(), len(off)))
    d0, d1, dd = mc.get("d_grid_mm", [-0.2, 0.2, 0.05])
    m0, m1, dm = mc.get("mf_grid", [0.4, 1.4, 0.1])
    dg = np.round(np.arange(d0, d1 + 1e-9, dd), 5)
    mg = np.round(np.arange(m0, m1 + 1e-9, dm), 5)
    print("channel tables: %d x %d 1D solutions ..." % (len(dg), len(mg)))
    T = MF.tables(case, dg, mg, C_hg=a.C_hg)
    P_h = case.Pin / 1e5
    nom = MF.network(T, np.zeros((1, len(off))), off, P_h)
    print("nominal channels: flow factor %.3f .. %.3f, max Twg %.0f C, outlet %.2f bar"
          % (nom["mf"].min(), nom["mf"].max(), nom["Twg"].max(), nom["Po"][0]))
    res = []
    for sig in mc.get("sigma_mm", [0.025, 0.05, 0.1]):
        for mu_ in mc.get("mean_mm", [0.0]):
            r = MF.monte_carlo(T, off, P_h, sig, mu_, mc.get("samples", 2000), mc.get("seed", 1))
            print("  MC sigma %.3f mean %+.3f: min flow p50 %.3f p01 %.3f | max Twg p50 %.0f p99 %.0f C | "
                  "min exit P p01 %.2f bar" % (sig, mu_, r["mf_min_p50"], r["mf_min_p01"], r["Twg_max_p50"],
                                               r["Twg_max_p99"], r["Pex_min_p01"]))
            res.append(r)
    plots.manifold(ang, off, nom["mf"][0], res, os.path.join(case.out_dir, "manifold.png"))
    _write_json(os.path.join(case.out_dir, "manifold.json"),
                dict(angles=ang, header_offset_bar=off, nominal_flow=nom["mf"][0], nominal_Twg=nom["Twg"][0],
                     monte_carlo=[{k: v for k, v in r.items() if not k.startswith("_")} for r in res]))


def cmd_radiation(case, a):
    from . import radiation as RD
    print(case.describe())
    res = case.couple()
    rc = case.cfg.get("radiation", {})
    if "X_H2O" in rc or "X_CO2" in rc:
        xw, xc = rc.get("X_H2O", 0.0), rc.get("X_CO2", 0.0)
    elif hasattr(case.gas, "table"):
        X = case.gas.table["chamber"]["X"]
        xw, xc = X.get("H2O", 0.0), X.get("CO2", 0.0)
    else:
        raise SystemExit("[radiation] needs X_H2O and X_CO2 with the ideal-gas model")
    ew = rc.get("eps_wall", 0.8)
    T_clip = rc.get("T_clip_K", 2500.0)
    GG = case.gas_grid()
    rads = {"clipped %.0f K" % T_clip: RD.estimate(res, GG, xw, xc, ew, T_clip)}
    if rc.get("extrapolate", True):
        rads["extrapolated"] = RD.estimate(res, GG, xw, xc, ew, None)
    print("radiation: X_H2O %.3f, X_CO2 %.3f, wall emissivity %.2f, beam length 0.9 D" % (xw, xc, ew))
    for lab, r in rads.items():
        print("  %-16s Q_rad %.2f kW (convective %.1f kW), max q_rad %.2f MW/m2, max q_rad/q %.3f, "
              "max dTwg %.1f K at x = %.1f mm -> max Twg %.0f C"
              % (lab, r["Q_rad_kW"], res["Qg_kW"], r["q_rad_max"] / 1e6, r["ratio_max"], r["dTwg_max"],
                 r["x_dTwg_max"], r["Twg_max_with_rad"]))
    plots.radiation(case, res, rads, os.path.join(case.out_dir, "radiation.png"))
    _write_json(os.path.join(case.out_dir, "radiation.json"),
                dict(X_H2O=xw, X_CO2=xc, eps_wall=ew, x=res["x"], cases=rads))


def cmd_injector(case, a):
    from . import injector as INJ
    print(case.describe())
    ic = case.cfg.get("injector")
    if not ic:
        raise SystemExit("no [injector] table in the case file")
    P2 = ic.get("P_chamber_bar", case.cfg["gas"]["pc_bar"]) * 1e5
    regen = None
    if any(isinstance(v, dict) and v.get("from_regen") for v in ic.values()):
        regen = case.regen_outlet(case.couple())
        print("regen outlet: T %.1f K, P %.2f bar" % (regen["Tout"], regen["P_out_bar"]))
    out = {}
    for name, e in ic.items():
        if not isinstance(e, dict):
            continue
        r = INJ.size(e, P2, case.fluid, regen)
        out[name] = r
        print("  %-10s %-8s T1 %.1f K, P1 %.2f bar -> %.2f bar: G %.0f kg/m2/s, Cd %.2f, %d x d %.3f mm%s"
              % (name, r["model"], r["T1_K"], r["P1_bar"], r["P2_bar"], r["G"], r["Cd"], r["n_holes"], r["d_mm"],
                 " (choked)" if r["info"].get("choked") else ""))
    _write_json(os.path.join(case.out_dir, "injector.json"), out)


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m regencool", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["gas", "run", "optimize", "tolerance", "manifold", "radiation", "injector"])
    p.add_argument("case", help="case TOML file")
    p.add_argument("-o", "--out", help="output folder (default results/<case name>)")
    p.add_argument("--C-hg", type=float, default=None, help="Bartz multiplier (overrides [case] C_hg)")
    p.add_argument("--no-stress", action="store_true", help="run: skip the 2D stations")
    p.add_argument("--axial", action="store_true", help="couple axial wall conduction (overrides [case])")
    a = p.parse_args(argv)
    case = Case.load(a.case, a.out)
    if a.C_hg is not None:
        case.C_hg = a.C_hg
    if a.axial:
        case.axial = True
    os.makedirs(case.out_dir, exist_ok=True)
    dict(gas=cmd_gas, run=cmd_run, optimize=cmd_optimize, tolerance=cmd_tolerance, manifold=cmd_manifold,
         radiation=cmd_radiation, injector=cmd_injector)[a.command](case, a)


if __name__ == "__main__":
    main()
