"""Atmosphere models. Each one provides density(h) in kg/m^3 for altitude h in m.

  Exponential(rho0, scale_height)  - rho0 * exp(-h / H)
  Table(altitudes, densities)      - smooth curve through log(density)
  us_standard_1976()               - Table built from the U.S. Standard Atmosphere 1976
  mars_pathfinder_1997()           - Table of the profile Pathfinder measured on entry
  titan_huygens_2005()             - Table of the profile Huygens measured on entry
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass(frozen=True)
class Exponential:
    rho0: float  # kg/m^3 at h = 0
    scale_height: float  # m

    def density(self, h):
        return self.rho0 * np.exp(-np.asarray(h) / self.scale_height)


class Table:
    """Density from a table: a cubic spline through log(density) vs altitude.

    A smooth curve (rather than straight segments) keeps the trajectory
    integrator fast. Outside the table it continues with the end slope.
    """

    def __init__(self, altitudes, densities, name: str = "table"):
        self.h = np.asarray(altitudes, dtype=float)
        self.log_rho = np.log(np.asarray(densities, dtype=float))
        if np.any(np.diff(self.h) <= 0):
            raise ValueError("altitudes must be strictly increasing")
        self.name = name
        self._spline = CubicSpline(self.h, self.log_rho)
        self._lo = float(self._spline(self.h[0], 1))
        self._hi = float(self._spline(self.h[-1], 1))

    def density(self, h):
        h = np.asarray(h, dtype=float)
        log_rho = self._spline(np.clip(h, self.h[0], self.h[-1]))
        log_rho = log_rho + self._lo * np.minimum(h - self.h[0], 0.0)
        log_rho = log_rho + self._hi * np.maximum(h - self.h[-1], 0.0)
        return np.exp(log_rho)

    def __repr__(self):
        return f"Table({self.name!r}, {self.h[0] / 1e3:.0f}-{self.h[-1] / 1e3:.0f} km)"


# --- U.S. Standard Atmosphere 1976 ------------------------------------------
# Below 86 km: the standard's own layer formulas (NOAA/NASA/USAF, 1976).
_R_AIR = 287.053  # J/kg/K
_G0 = 9.80665  # m/s^2
_R0 = 6356.766e3  # m, radius used for geopotential altitude
_LAYERS = [  # geopotential base (m), lapse rate (K/m)
    (0.0, -6.5e-3), (11e3, 0.0), (20e3, 1.0e-3), (32e3, 2.8e-3),
    (47e3, 0.0), (51e3, -2.8e-3), (71e3, -2.0e-3),
]
# 86 km and above: densities from the standard's published tables.
_UPPER = [  # geometric altitude (km), density (kg/m^3)
    (86, 6.958e-6), (90, 3.416e-6), (95, 1.393e-6), (100, 5.604e-7),
    (105, 2.325e-7), (110, 9.708e-8), (115, 4.289e-8), (120, 2.222e-8),
    (125, 1.291e-8), (130, 8.152e-9), (135, 5.465e-9), (140, 3.831e-9),
    (145, 2.781e-9), (150, 2.076e-9),
]


def _us76_lower(z):
    """Density below 86 km at geometric altitude z (m)."""
    hgp = _R0 * z / (_R0 + z)
    t, p = 288.15, 101325.0
    for i, (base, lapse) in enumerate(_LAYERS):
        top = _LAYERS[i + 1][0] if i + 1 < len(_LAYERS) else np.inf
        dh = min(hgp, top) - base
        if lapse == 0.0:
            p *= np.exp(-_G0 * dh / (_R_AIR * t))
        else:
            t_new = t + lapse * dh
            p *= (t / t_new) ** (_G0 / (_R_AIR * lapse))
            t = t_new
        if hgp <= top:
            break
    return p / (_R_AIR * t)


def us_standard_1976(step: float = 1000.0) -> Table:
    low = np.arange(0.0, 86e3, step)
    h = np.concatenate([low, [a * 1e3 for a, _ in _UPPER]])
    rho = np.concatenate([[_us76_lower(z) for z in low], [r for _, r in _UPPER]])
    return Table(h, rho, name="US Standard 1976")


def mars_pathfinder_1997() -> Table:
    """Mars atmosphere as measured by Pathfinder during its entry, 4 July 1997.

    NASA PDS data set MPFL-M-ASIMET-4-DDR-EDL-V1.0 (see the data file header).
    Covers 9 to 141 km above the landing site, whose radius is 3389.72 km.
    The densities were derived from the capsule's own deceleration, so
    matching Pathfinder's deceleration with this table is a consistency
    check, not an independent test.
    """
    path = Path(__file__).parent / "data" / "mars_pathfinder_atmosphere.csv"
    alt_km, rho, _, _ = np.loadtxt(path, delimiter=",", skiprows=6).T
    return Table(alt_km[::-1] * 1e3, rho[::-1], name="Mars Pathfinder 1997")


def titan_huygens_2005(molar_mass: float = 28.0) -> Table:
    """Titan atmosphere as measured by Huygens (HASI), 14 January 2005.

    NASA PDS / ESA PSA data set HP-SSA-HASI-2-3-4-MISSION-V1.1 (see the data
    file header). Covers the surface to 1380 km. Below 147 km the archive
    gives density directly. Above 157 km it gives only pressure and
    temperature, so density is computed as p M / (R T) with `molar_mass` in
    g/mol; the archived descent rows imply 27.8 to 28.6, so expect about 3%
    uncertainty from this step. As with Pathfinder, the upper profile was
    derived from the probe's own deceleration.
    """
    path = Path(__file__).parent / "data" / "titan_huygens_atmosphere.csv"
    alt, rho = [], []
    for line in path.read_text().splitlines():
        if not line[:1].isdigit():
            continue
        a, d, p, t = line.split(",")
        alt.append(float(a) * 1e3)
        rho.append(float(d) if d else float(p) * molar_mass / (8314.46 * float(t)))
    return Table(alt[::-1], rho[::-1], name="Titan Huygens 2005")
