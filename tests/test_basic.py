import math
import os

import numpy as np
import pytest

from regencool.case import Case
from regencool.contour import Contour
from regencool.fe2d import Material, mesh, structural
from regencool.gas import mach_from_area

HERE = os.path.dirname(os.path.abspath(__file__))
DEMO = os.path.join(HERE, "..", "examples", "demo_engine.toml")


def test_parametric_contour():
    c = Contour.parametric(22.0, 9.0, 70.0, eps_exit=5.0)
    assert abs(c.r_t - 9.0) < 1e-3
    assert abs(c.eps_exit - 5.0) < 1e-3
    assert c.x0 == 0.0 and c.x0 < c.x_t < c.x_end
    assert abs(float(c.R(10.0)) - 22.0) < 1e-9


def test_mach_area_round_trip():
    g = 1.22
    for M in (0.2, 0.6, 1.5, 3.0):
        ar = (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1)))
        assert abs(float(mach_from_area(ar, g, M > 1)[0]) - M) < 1e-6


@pytest.fixture(scope="module")
def demo_case(tmp_path_factory):
    return Case.load(DEMO, out_dir=str(tmp_path_factory.mktemp("demo")), log=lambda *a: None)


def test_energy_balance(demo_case):
    res = demo_case.solve()
    # heat through the hot wall = coolant enthalpy rise
    assert abs(res["Q_kW"] * 1e3 - demo_case.mdot * res["dh_coolant"]) / (res["Q_kW"] * 1e3) < 0.01
    m = demo_case.metrics(res)
    assert m["dP_ch"] > 0 and m["dP"] > m["dP_ch"]
    assert m["Tout"] > demo_case.Tin - 273.15


def test_free_thermal_expansion_is_stress_free():
    """Uniform temperature rise, no channel, no pressure: a sector of a free ring expands without stress."""
    mat = Material()
    M = mesh(10.0, 13.0, 2 * math.pi / 24, 1.0, "tear", 1.0, h=0.05, hmax=0.1, channel=False)
    T = np.full(len(M["XN"]), mat.T_ref + 500.0)
    res, _ = structural(M, T, 0.0, 0.0, 0.0, mat)
    scale = mat.E * mat.alpha * 500.0
    assert res["tot"]["vm"].max() < 1e-3 * scale
    assert abs(res["tot"]["ez"] - mat.alpha * 500.0) < 1e-6


def test_cli_run(tmp_path):
    from regencool.__main__ import main
    main(["run", DEMO, "-o", str(tmp_path), "--no-stress"])
    assert (tmp_path / "summary.json").exists()
    assert (tmp_path / "solution_1d.csv").exists()


def test_cantera_chamber():
    pytest.importorskip("cantera")
    from regencool import cantera_gas
    fuel = dict(formula="C3H8O", hf_kJ_mol=-318.1, cp_J_kgK=2600.0, T_K=293.15)
    ox = dict(species="N2O", T_K=293.15, dh_vap_kJ_kg=170.0)
    t = cantera_gas.compute(25e5, 4.0, fuel, ox, 6.0, npts=200, modes=("eq",))
    assert 2800 < t["chamber"]["T"] < 3300
    assert 1400 < t["cstar_eq"] < 1700
    eps = [row["eps"] for row in t["exp_eq"]]
    assert max(eps) >= 6.0


def test_axial_conduction(demo_case):
    demo_case.axial = True
    try:
        res = demo_case.solve()
    finally:
        demo_case.axial = False
    assert res["axial"] and res["axial_conv"]
    # axial heat only moves between stations: gas-side total = coolant pickup
    assert abs(res["Qg_kW"] - res["Q_kW"]) / res["Q_kW"] < 0.02
    assert abs(res["Q_kW"] * 1e3 - demo_case.mdot * res["dh_coolant"]) / (res["Q_kW"] * 1e3) < 0.01
    # spreading can only lower the peak
    assert res["Twg"].max() <= res["Twg_noax"].max() + 1.0


def test_ipa_saturation():
    from regencool import boiling
    sat = boiling.ipa_saturation()
    s = sat.at(101325.0)
    assert abs(s["Tsat"] - 355.4) < 1.0
    assert 600e3 < s["hfg"] < 720e3
    assert s["rho_l"] > 100 * s["rho_v"]
    assert sat.at(20e5)["Tsat"] > s["Tsat"]


def test_boiling_check(demo_case):
    from regencool import boiling
    res = demo_case.solve()
    assert demo_case.boiling(res) is None            # demo channels stay supercritical
    sat = boiling.ipa_saturation()
    hot = dict(res, P=np.full_like(res["P"], 10e5))  # same temperatures at 10 bar -> wall far above Tsat
    b = boiling.check(hot, sat)
    assert np.isfinite(b["Tsat_margin"]).all() and b["flag"].max() == 2
    cold = dict(hot, Twc=np.minimum(res["Tb"], 330.0), Tb=np.minimum(res["Tb"], 330.0))
    assert boiling.check(cold, sat)["flag"].max() == 0


def test_injector_models():
    from regencool import injector
    G, info = injector.ipa_real(293.15, 30e5, 20e5)
    G_spi = injector.spi(info["rho1"], 30e5, 20e5)
    assert abs(G / G_spi - 1) < 0.02 and not info["choked"]
    G_hot, info_hot = injector.ipa_real(495.0, 55e5, 25e5)
    assert info_hot["choked"] and G_hot < injector.spi(info_hot["rho1"], 55e5, 25e5)
    Gn, info_n = injector.n2o_dyer(273.15, 36e5, 25e5)
    assert info_n["G_spi"] > Gn > 0
    assert abs(injector.diameter(1.0, 1e4, 1.0, 1) - math.sqrt(4e-4 / math.pi) * 1e3) < 1e-9


def test_radiation_estimate(demo_case):
    from regencool import radiation
    res = demo_case.solve()
    r = radiation.estimate(res, demo_case.gas_grid(), 0.3, 0.15)
    assert 0 < r["Q_rad_kW"] < 0.5 * res["Qg_kW"]
    assert (r["eps_g"] > 0).all() and (r["eps_g"] < 1).all()
    assert r["dTwg_max"] > 0


def test_cli_radiation_injector(tmp_path):
    from regencool.__main__ import main
    main(["radiation", DEMO, "-o", str(tmp_path)])
    main(["injector", DEMO, "-o", str(tmp_path)])
    assert (tmp_path / "radiation.json").exists()
    assert (tmp_path / "injector.json").exists()
