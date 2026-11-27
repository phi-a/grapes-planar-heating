"""Monte Carlo dispersion: how far the outputs move when the inputs vary.

Four inputs are perturbed, each drawn independently from a normal
distribution. The default 1-sigma values are assumptions, typical of entry
dispersion analyses, not mission-specific figures.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .planar import simulate

SIGMA = dict(
    gamma=np.radians(0.05),  # entry flight-path angle, rad
    speed=0.001,  # entry speed, fraction
    drag=0.02,  # drag coefficient, fraction
    density=0.05,  # atmospheric density at every altitude, fraction
)


@dataclass(frozen=True)
class _Scaled:
    """An atmosphere with its density multiplied by a constant."""

    base: object
    factor: float

    def density(self, h):
        return self.factor * self.base.density(h)


def disperse(entry, t_end, n=400, seed=0, sigma=SIGMA):
    """Fly `n` perturbed copies of `entry`.

    Returns arrays of peak deceleration (m/s^2) and altitude at t_end (m),
    with the unperturbed run first.
    """
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((n, 4))
    z[0] = 0.0
    peak, alt = np.empty(n), np.empty(n)
    for i, (zg, zv, zd, zr) in enumerate(z):
        body = entry["body"]
        run = {
            **entry,
            "gamma0": entry["gamma0"] + sigma["gamma"] * zg,
            "v0": entry["v0"] * (1 + sigma["speed"] * zv),
            "vehicle": replace(entry["vehicle"], cd=entry["vehicle"].cd * (1 + sigma["drag"] * zd)),
            "body": replace(body, atmosphere=_Scaled(body.atmosphere, 1 + sigma["density"] * zr)),
        }
        tr = simulate(**run, t_end=t_end, dt_out=0.1, rtol=1e-8)
        peak[i], alt[i] = tr.decel.max(), tr.h[-1]
    return peak, alt
