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
