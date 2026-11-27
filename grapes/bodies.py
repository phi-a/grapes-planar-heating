"""Planetary bodies: size, gravity, rotation and an atmosphere model."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from .atmosphere import Exponential, titan_huygens_2005, us_standard_1976
from .heating import (TAUBER_SUTTON_EARTH, TAUBER_SUTTON_MARS,
                      TITAN_HUYGENS_NONBOLTZMANN, sutton_graves_coefficient)


@dataclass(frozen=True)
class Body:
    name: str
    radius: float  # m, used as the reference sphere
    mu: float  # m^3/s^2, gravitational parameter
    atmosphere: Any  # anything with .density(h); see atmosphere.py
    k_sutton_graves: float  # sqrt(kg)/m, convective coefficient for this gas mix
    rotation_rate: float = 0.0  # rad/s; used by frames.py and simulate(latitude=, azimuth=)
    radiative: Any = None  # anything with .rate(rho, v, nose_radius); see heating.py

    def density(self, h):
        return self.atmosphere.density(h)

    def gravity(self, h):
        return self.mu / (self.radius + np.asarray(h)) ** 2

    def with_atmosphere(self, atmosphere) -> "Body":
        """Same body, different atmosphere model."""
        return replace(self, atmosphere=atmosphere)


# Convective coefficients from Sutton & Graves (1971); see heating.py.
EARTH = Body(
    name="Earth",
    radius=6378.1e3,
    mu=3.9860e14,
    atmosphere=us_standard_1976(),
    k_sutton_graves=sutton_graves_coefficient(dict(N2=0.79, O2=0.21)),
    rotation_rate=7.2921066e-5,
    radiative=TAUBER_SUTTON_EARTH,
)

MARS = Body(
    name="Mars",
    radius=3396.2e3,
    mu=4.2830e13,
    atmosphere=Exponential(rho0=0.02, scale_height=11100.0),
    k_sutton_graves=sutton_graves_coefficient(dict(CO2=0.97, N2=0.03)),
    rotation_rate=7.0882360e-5,
    radiative=TAUBER_SUTTON_MARS,
)

# Titan. The radiative model is a three-point fit for the Huygens probe only
# (see heating.PowerLawFit); swap it with dataclasses.replace for the
# Boltzmann upper estimate.
# k_sutton_graves: from Sutton & Graves (1971) for 96.99% N2, 2.30% CH4,
# 0.71% Ar by mole, the composition the Huygens heating study used.
TITAN = Body(
    name="Titan",
    radius=2575.0e3,
    mu=8.9781e12,
    atmosphere=titan_huygens_2005(),
    k_sutton_graves=sutton_graves_coefficient(dict(N2=0.9699, CH4=0.0230, Ar=0.0071)),
    rotation_rate=4.5607e-6,  # 15.945-day rotation
    radiative=TITAN_HUYGENS_NONBOLTZMANN,
)


BODIES = {b.name.lower(): b for b in (EARTH, MARS, TITAN)}
