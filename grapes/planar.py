"""Planar (2-D) entry over a spherical body.

State: altitude h, speed V, flight-path angle gamma, range angle theta,
integrated convective heat load Q. Speed and angle are relative to the
atmosphere.

Planet rotation is optional. Give simulate() the entry latitude and heading
and it adds the Coriolis and centrifugal terms that act in the entry plane,
holding latitude and heading at their entry values. Leave them out for a
non-rotating planet. Turning out of the plane needs the 3-D model.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .bodies import Body
from .heating import CONVECTIVE
from .vehicle import Vehicle


@dataclass
class Trajectory:
    t: np.ndarray  # s
    h: np.ndarray  # m
    v: np.ndarray  # m/s
    gamma: np.ndarray  # rad
    range_angle: np.ndarray  # rad
    heat_load: np.ndarray  # J/m^2, convective
    q_conv: np.ndarray  # W/m^2
    q_rad: np.ndarray  # W/m^2; zero with no radiative model, NaN where the model has no data
    decel: np.ndarray  # m/s^2, drag only
    body: Body
    vehicle: Vehicle
    skipped: bool = False  # bounced back out above the entry altitude

    @property
    def range(self) -> np.ndarray:
        """Ground range along the reference sphere, m."""
        return self.body.radius * self.range_angle

    @property
    def q_total(self) -> np.ndarray:
        """Convective plus radiative, counting radiative as zero where it has no data."""
        return self.q_conv + np.nan_to_num(self.q_rad)

    @property
    def heat_load_rad(self) -> np.ndarray:
        """Integrated radiative heat load, J/m^2."""
        q = np.nan_to_num(self.q_rad)
        steps = 0.5 * (q[1:] + q[:-1]) * np.diff(self.t)
        return np.concatenate([[0.0], np.cumsum(steps)])

    @property
    def peak_q_conv(self) -> float:
        return float(self.q_conv.max())

    @property
    def peak_q_rad(self) -> float:
        return float(np.nan_to_num(self.q_rad).max())

    @property
    def peak_q_total(self) -> float:
        return float(self.q_total.max())

    @property
    def peak_decel_g(self) -> float:
        return float(self.decel.max() / 9.80665)


def _rates(t, x, body, veh, heat, rot):
    h, v, gam, _theta, _q = x
    rho = body.density(h)
    g = body.gravity(h)
    r = body.radius + h
    bc = veh.ballistic_coefficient
    drag_acc = rho * v**2 / (2.0 * bc)
    sg, cg = np.sin(gam), np.cos(gam)
    v_dot = -drag_acc - g * sg
    gam_dot = drag_acc * veh.lift_to_drag / v - (g / v - v / r) * cg
    if rot is not None:
        # Rotating-planet terms (e.g. Vinh, Busemann & Culp, 1980), with
        # latitude and heading (clockwise from north) frozen at entry.
        w, clat, slat, spsi, cpsi = rot
        v_dot += w**2 * r * clat * (sg * clat - cg * slat * cpsi)
        gam_dot += 2.0 * w * clat * spsi + w**2 * r / v * clat * (cg * clat + sg * slat * cpsi)
    return [
        v * sg,
        v_dot,
        gam_dot,
        v * np.cos(gam) / r,
        heat(rho, v, veh.nose_radius, body),
    ]


def simulate(
    body: Body,
    vehicle: Vehicle,
    h0: float,
    v0: float,
    gamma0: float,
    t_end: float,
    dt_out: float = 1.0,
    convective: str = "sutton_graves",
    stop_altitude: float = 0.0,
    rtol: float = 1e-10,
    atol: float = 1e-8,
    latitude: float | None = None,
    azimuth: float | None = None,
) -> Trajectory:
    """Fly one entry. gamma0 in radians; v0 and gamma0 relative to the atmosphere.

    Give `latitude` and `azimuth` (radians, heading clockwise from north,
    relative to the atmosphere) to include planet rotation; see the module
    docstring.

    Stops at t_end, at stop_altitude, or when the vehicle skips back out
    above its entry altitude.
    """
    heat = CONVECTIVE[convective]
    if (latitude is None) != (azimuth is None):
        raise ValueError("give both latitude and azimuth, or neither")
    rot = None
    if latitude is not None:
        rot = (body.rotation_rate, np.cos(latitude), np.sin(latitude),
               np.sin(azimuth), np.cos(azimuth))

    def skip_out(t, x, *args):
        return x[0] - (h0 + 1.0)

    skip_out.terminal = True
    skip_out.direction = 1

    def hit_floor(t, x, *args):
        return x[0] - stop_altitude

    hit_floor.terminal = True
    hit_floor.direction = -1

    t_eval = np.arange(0.0, t_end + 0.5 * dt_out, dt_out)
    sol = solve_ivp(
        _rates,
        (0.0, t_eval[-1]),
        [h0, v0, gamma0, 0.0, 0.0],
        method="DOP853",
        t_eval=t_eval,
        events=(hit_floor, skip_out),
        args=(body, vehicle, heat, rot),
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    h, v, gam, theta, q_load = sol.y
    rho = body.density(h)
    if body.radiative is None:
        q_rad = np.zeros_like(h)
    else:
        q_rad = body.radiative.rate(rho, v, vehicle.nose_radius)
    return Trajectory(
        t=sol.t,
        h=h,
        v=v,
        gamma=gam,
        range_angle=theta,
        heat_load=q_load,
        q_conv=heat(rho, v, vehicle.nose_radius, body),
        q_rad=q_rad,
        decel=rho * v**2 / (2.0 * vehicle.ballistic_coefficient),
        body=body,
        vehicle=vehicle,
        skipped=len(sol.t_events[1]) > 0,
    )
