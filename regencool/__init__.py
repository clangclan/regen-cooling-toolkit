"""regencool: 1D regenerative-cooling and 2D channel-section thermo-elastic toolkit for small liquid rocket chambers.

Main entry points
-----------------
Case            -- loads a TOML case file and builds contour, gas model, coolant and channel design
solver1d.solve  -- 1D coupled hot-gas / wall / coolant march
station2d.station -- 2D thermo-elastic analysis of one periodic channel sector at an axial station
python -m regencool {gas,run,optimize,tolerance,manifold} case.toml
"""
__version__ = "0.1.0"

from .case import Case  # noqa: E402,F401
