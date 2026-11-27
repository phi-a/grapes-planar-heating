"""Sweep entry speed and flight-path angle to map the entry corridor."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bodies import Body
from .planar import simulate
from .vehicle import Vehicle


@dataclass
class CorridorMap:
    speeds: np.ndarray  # m/s, shape (n_v,)
    angles: np.ndarray  # rad, shape (n_g,)
    peak_q: np.ndarray  # W/m^2, convective, shape (n_g, n_v)
    peak_q_rad: np.ndarray  # W/m^2
    peak_q_total: np.ndarray  # W/m^2, peak of convective + radiative
    peak_g: np.ndarray  # Earth g
    heat_load: np.ndarray  # J/m^2
    skipped: np.ndarray  # bool


def entry_corridor(
    body: Body,
    vehicle: Vehicle,
    h0: float,
    speeds,
    angles,
    convective: str = "sutton_graves",
    t_end: float = 1500.0,
    stop_altitude: float = 20e3,
    convert=None,
    **sim_kw,
) -> CorridorMap:
    """Fly every (speed, angle) pair. Angles in radians, negative = descending.

    `convert(v, gamma) -> (v, gamma)` is applied to each pair before flying,
    e.g. to sweep space-fixed entry states (see frames.inertial_to_relative).
    Extra keywords (latitude, azimuth) go to simulate().
    """
    speeds = np.asarray(speeds, dtype=float)
    angles = np.asarray(angles, dtype=float)
    shape = (angles.size, speeds.size)
    peak_q, peak_rad, peak_tot, peak_g, load = (np.empty(shape) for _ in range(5))
    skipped = np.zeros(shape, dtype=bool)
    for i, gam in enumerate(angles):
        for j, v in enumerate(speeds):
            vr, gr = convert(v, gam) if convert else (v, gam)
            tr = simulate(
                body, vehicle, h0, vr, gr, t_end,
                dt_out=0.25, convective=convective,
                stop_altitude=stop_altitude, rtol=1e-8, **sim_kw,
            )
            peak_q[i, j] = tr.peak_q_conv
            peak_rad[i, j] = tr.peak_q_rad
            peak_tot[i, j] = tr.peak_q_total
            peak_g[i, j] = tr.peak_decel_g
            load[i, j] = tr.heat_load[-1]
            skipped[i, j] = tr.skipped
    return CorridorMap(speeds, angles, peak_q, peak_rad, peak_tot, peak_g, load, skipped)


def skip_out_boundary(body, vehicle, h0, speeds, shallow, steep, tol=1e-4,
                      convert=None, **kw):
    """Shallowest angle (rad) that still captures, per speed, by bisection.

    Returns NaN where even the `shallow` angle captures.
    """
    def skips(v, gam):
        if convert:
            v, gam = convert(v, gam)
        return simulate(body, vehicle, h0, v, gam, 1500.0, dt_out=1.0,
                        stop_altitude=20e3, rtol=1e-8, **kw).skipped

    out = np.full(len(speeds), np.nan)
    for j, v in enumerate(speeds):
        if not skips(v, shallow):
            continue
        a, b = shallow, steep  # a skips, b captures
        while abs(a - b) > tol:
            mid = 0.5 * (a + b)
            a, b = (mid, b) if skips(v, mid) else (a, mid)
        out[j] = b
    return out
