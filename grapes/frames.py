"""Entry-state conversions between the space-fixed and rotating frames.

The planar model flies relative to the atmosphere and knows nothing about
planet rotation. Published entry states are usually space-fixed ("inertial"),
so convert once, up front:

    v, gamma, _ = inertial_to_relative(EARTH, h0, v_inertial, gamma_inertial,
                                       latitude, azimuth)
    simulate(EARTH, vehicle, h0, v, gamma, ...)

This captures the main effect of rotation on heating and loads (the airspeed
at entry). It does not capture Coriolis drift during the entry; that needs
the 3-D model.
"""
from __future__ import annotations

import numpy as np

from .bodies import Body


def inertial_to_relative(body: Body, h, v, gamma, latitude, azimuth):
    """Space-fixed entry state -> state relative to the rotating atmosphere.

    Angles in radians. Azimuth is measured clockwise from north.
    Returns (speed, flight-path angle, azimuth) relative to the atmosphere.
    """
    east = v * np.cos(gamma) * np.sin(azimuth)
    north = v * np.cos(gamma) * np.cos(azimuth)
    up = v * np.sin(gamma)
    east = east - body.rotation_rate * (body.radius + h) * np.cos(latitude)
    v_rel = np.sqrt(east**2 + north**2 + up**2)
    return v_rel, np.arcsin(up / v_rel), np.arctan2(east, north) % (2 * np.pi)
