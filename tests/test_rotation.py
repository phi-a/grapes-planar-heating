"""Rotating-planet terms in the planar model."""
from dataclasses import replace

import numpy as np
import pytest

from grapes import EARTH, simulate
from grapes.cases import STARDUST

ENTRY = dict(vehicle=STARDUST, h0=125e3, v0=11e3, gamma0=np.radians(-7.0), t_end=120)


def test_zero_spin_equals_no_rotation():
    still = replace(EARTH, rotation_rate=0.0)
    a = simulate(still, **ENTRY)
    b = simulate(still, **ENTRY, latitude=0.3, azimuth=1.0)
    np.testing.assert_allclose(a.h, b.h, rtol=1e-12)
    np.testing.assert_allclose(a.v, b.v, rtol=1e-12)


def test_eastward_is_lifted_and_westward_is_pushed_down():
    none = simulate(EARTH, **ENTRY)
    east = simulate(EARTH, **ENTRY, latitude=0.0, azimuth=np.pi / 2)
    west = simulate(EARTH, **ENTRY, latitude=0.0, azimuth=3 * np.pi / 2)
    assert east.peak_decel_g < none.peak_decel_g < west.peak_decel_g


def test_latitude_and_azimuth_must_come_together():
    with pytest.raises(ValueError):
        simulate(EARTH, **ENTRY, latitude=0.1)
