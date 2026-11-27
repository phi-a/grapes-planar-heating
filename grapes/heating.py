"""Stagnation-point heating models. All rates in W/m^2, SI inputs.

Convective (both the same form, q = k * sqrt(rho / Rn) * V^3):
  - chapman:        Earth-only classic, k = 1.83e-4. Kept for comparison.
  - sutton_graves:  k depends on the atmosphere's gas mix (see Body).

Radiative:
  - TauberSutton:   q = C * Rn^a * rho^b * f(V), one coefficient set per
                    atmosphere. Earth and Mars sets are below.
  - PowerLawFit:    q = C * rho^a * V^b through published reference points,
                    for bodies with no general fit (Titan).
Any object with a .rate(rho, v, nose_radius) method works as a radiative model.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence, Union

import numpy as np

K_CHAPMAN = 1.83e-4  # sqrt(kg)/m


def _convective(k, rho, v, nose_radius):
    return k * np.sqrt(np.asarray(rho) / nose_radius) * np.asarray(v) ** 3


def chapman(rho, v, nose_radius, body=None):
    return _convective(K_CHAPMAN, rho, v, nose_radius)


def sutton_graves(rho, v, nose_radius, body):
    return _convective(body.k_sutton_graves, rho, v, nose_radius)


CONVECTIVE = {"chapman": chapman, "sutton_graves": sutton_graves}


# Sutton & Graves, "A General Stagnation-Point Convective-Heating Equation for
# Arbitrary Gas Mixtures", NASA TR R-376, 1971 (NTRS 19720003329).
# Heat-transfer coefficient K of each base gas, kg/(s m^1.5 atm^0.5) (Table II),
# and molecular weight (Table I).
SUTTON_GRAVES_GASES = {
    "N2": (0.1112, 28.014), "O2": (0.1201, 32.000), "H2": (0.0395, 2.016),
    "He": (0.0797, 4.003), "Ne": (0.1474, 20.183), "Ar": (0.1495, 39.948),
    "CO2": (0.1210, 44.011), "NH3": (0.0990, 17.031), "CH4": (0.0807, 16.043),
}


def sutton_graves_coefficient(mole_fractions: dict) -> float:
    """k for q = k sqrt(rho / Rn) V^3, in sqrt(kg)/m, for a gas mixture.

    Mixes the base-gas coefficients with the report's eq. 44,
    1/K^2 = sum(c_i / K_i^2) over mass fractions c_i, then converts K to k
    using a stagnation pressure of rho V^2 and an enthalpy of V^2 / 2:
    k = K / (2 sqrt(101325)).
    """
    mass = {g: x * SUTTON_GRAVES_GASES[g][1] for g, x in mole_fractions.items()}
    total = sum(mass.values())
    inv_k2 = sum(m / total / SUTTON_GRAVES_GASES[g][0] ** 2 for g, m in mass.items())
    return float(inv_k2 ** -0.5 / (2.0 * np.sqrt(101325.0)))


@dataclass(frozen=True)
class TauberSutton:
    """Tauber & Sutton radiative fit: q = C * Rn^a * rho^b * f(V), in W/cm^2
    with Rn in m, rho in kg/m^3, V in m/s. `rate` returns W/m^2.

    f(V) is a table, interpolated linearly as the paper says. Below the
    table's lowest speed the rate is zero; above its highest speed f is held
    at the last value, which under-predicts, so stay inside `v_range`.

    `valid` records the paper's stated ranges (speed m/s, density kg/m^3,
    nose radius m). Nothing enforces them.
    """

    name: str
    c: float
    a: Union[float, Callable]  # constant, or a(rho, v, nose_radius)
    b: float
    f_speeds: Sequence[float]  # m/s
    f_values: Sequence[float]
    valid: dict = None

    @property
    def v_range(self):
        return self.f_speeds[0], self.f_speeds[-1]

    def rate(self, rho, v, nose_radius):
        rho, v = np.asarray(rho, dtype=float), np.asarray(v, dtype=float)
        f = np.interp(v, self.f_speeds, self.f_values)
        f = np.where(v < self.f_speeds[0], 0.0, f)
        vv = np.clip(v, *self.v_range)
        a = self.a(rho, vv, nose_radius) if callable(self.a) else self.a
        return 1e4 * self.c * nose_radius**a * rho**self.b * f


def _earth_exponent(rho, v, nose_radius):
    a = 1.072e6 * v**-1.88 * rho**-0.325
    cap = 1.0 if nose_radius < 1.0 else 0.6 if nose_radius <= 2.0 else 0.5
    return np.minimum(a, cap)


# Tauber & Sutton, "Stagnation-Point Radiative Heating Relations for Earth and
# Mars Entries", J. Spacecraft and Rockets 28(1), 1991, pp. 40-42.
# Eq. 2 (Earth), Eq. 3 (Mars) and Table 1 (f); checked against the paper.
# Cold wall, no ablation, thermochemical equilibrium. The paper reports
# agreement with detailed computations within about 20% over the valid range.
TAUBER_SUTTON_EARTH = TauberSutton(
    name="Tauber-Sutton, Earth (air)",
    c=4.736e4,
    a=_earth_exponent,
    b=1.22,
    f_speeds=[9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750, 11000, 11500,
              12000, 12500, 13000, 13500, 14000, 14500, 15000, 15500, 16000],
    f_values=[1.5, 4.3, 9.7, 19.5, 35, 55, 81, 115, 151, 238,
              359, 495, 660, 850, 1065, 1313, 1550, 1780, 2040],
    valid=dict(speed=(10e3, 16e3), density=(6.66e-5, 6.31e-4), nose_radius=(0.3, 3.0)),
)

TAUBER_SUTTON_MARS = TauberSutton(
    name="Tauber-Sutton, Mars (CO2-N2)",
    c=2.35e4,
    a=0.526,
    b=1.19,
    f_speeds=[6000, 6150, 6300, 6500, 6700, 6900, 7000, 7200, 7400, 7600,
              7800, 8000, 8200, 8400, 8600, 8800, 9000],
    f_values=[0.2, 1.0, 1.95, 3.42, 5.1, 7.1, 8.1, 10.2, 12.5, 14.8,
              17.1, 19.2, 21.4, 24.1, 26.0, 28.9, 32.8],
    valid=dict(speed=(6.5e3, 9e3), density=(1e-4, 1e-3), nose_radius=(1.0, 23.0)),
)


class PowerLawFit:
    """Radiative fit q = C * rho^a * V^b through three published points.

    Three points fix the three constants exactly, so this reproduces the
    source at those points and is only an interpolation between them.
    Outside the speed span of the points it returns NaN ("no data") rather
    than extrapolate. It has no nose-radius dependence: use it only for the
    vehicle the points were computed for.
    """

    def __init__(self, name, points, nose_radius):
        """points: three (speed m/s, density kg/m^3, q W/cm^2)."""
        self.name, self.points, self.nose_radius = name, tuple(points), nose_radius
        v, rho, q = (np.array(c, dtype=float) for c in zip(*self.points))
        a_mat = np.column_stack([np.ones(3), np.log(rho), np.log(v)])
        self.log_c, self.a, self.b = np.linalg.solve(a_mat, np.log(q))
        self.valid = dict(speed=(v.min(), v.max()), density=(rho.min(), rho.max()),
                          nose_radius=(nose_radius, nose_radius))

    def rate(self, rho, v, nose_radius=None):
        rho, v = np.asarray(rho, dtype=float), np.asarray(v, dtype=float)
        q = 1e4 * np.exp(self.log_c + self.a * np.log(rho) + self.b * np.log(v))
        lo, hi = self.valid["speed"]
        return np.where((v >= lo) & (v <= hi), q, np.nan)


# Huygens at Titan: stagnation-point values near peak heating from Johnston,
# Hollis & Sutton, "Radiative Heating Methodology for the Huygens Probe",
# AIAA 2006-3426, Tables 1, 4, 5 and 8. Times are seconds on the pre-flight
# design trajectory. Heat fluxes in W/cm^2, flowfield-radiation coupled.
HUYGENS_DESIGN_POINTS = dict(
    t=(185.0, 189.0, 193.0),
    v=(5489.0, 5126.0, 4705.0),  # m/s
    rho=(1.83e-4, 2.96e-4, 3.79e-4),  # kg/m^3
    q_conv=(37.6, 38.7, 33.7),
    q_conv_uncoupled=(43.9, 44.5, 37.9),
    q_rad_boltzmann=(59.0, 66.5, 54.0),  # CN states in equilibrium: upper estimate
    q_rad_nonboltzmann=(20.8, 32.8, 34.9),  # collisional-radiative CN model: lower estimate
)
_P = HUYGENS_DESIGN_POINTS

TITAN_HUYGENS_NONBOLTZMANN = PowerLawFit(
    "Titan, Huygens points (non-Boltzmann CN)",
    zip(_P["v"], _P["rho"], _P["q_rad_nonboltzmann"]), nose_radius=1.25)
TITAN_HUYGENS_BOLTZMANN = PowerLawFit(
    "Titan, Huygens points (Boltzmann CN)",
    zip(_P["v"], _P["rho"], _P["q_rad_boltzmann"]), nose_radius=1.25)
