import numpy as np
import pytest

from grapes import EARTH, simulate
from grapes.atmosphere import Exponential
from grapes.bodies import Body
from grapes.cases import STARDUST_ENTRY, STARDUST_FLIGHT
from grapes.frames import inertial_to_relative

H, V, GAM = 120e3, 11e3, np.radians(-6.0)


def test_no_rotation_changes_nothing():
    still = Body("still", EARTH.radius, EARTH.mu, Exponential(1.2, 8500.0), 1.7e-4)
    v, gam, az = inertial_to_relative(still, H, V, GAM, 0.3, 1.0)
    assert (v, gam, az) == pytest.approx((V, GAM, 1.0))


def test_due_east_on_equator_loses_ground_speed():
    v, gam, az = inertial_to_relative(EARTH, H, V, 0.0, 0.0, np.pi / 2)
    assert v == pytest.approx(V - EARTH.rotation_rate * (EARTH.radius + H))
    assert gam == pytest.approx(0.0)
    assert az == pytest.approx(np.pi / 2)


def test_due_west_gains_it():
    v, _, _ = inertial_to_relative(EARTH, H, V, 0.0, 0.0, 3 * np.pi / 2)
    assert v == pytest.approx(V + EARTH.rotation_rate * (EARTH.radius + H))


def test_pole_is_unaffected():
    v, gam, _ = inertial_to_relative(EARTH, H, V, GAM, np.pi / 2, 1.0)
    assert (v, gam) == pytest.approx((V, GAM))


def test_vertical_speed_is_preserved():
    v, gam, _ = inertial_to_relative(EARTH, H, V, GAM, 0.7, 1.8)
    assert v * np.sin(gam) == pytest.approx(V * np.sin(GAM))


def test_stardust_peak_g_near_flight():
    # 3.3% low: standard (not that night's) atmosphere, constant drag.
    tr = simulate(**STARDUST_ENTRY, t_end=140, dt_out=0.1)
    assert tr.peak_decel_g == pytest.approx(STARDUST_FLIGHT["peak_g"], rel=0.06)


def test_stardust_altitude_at_drogue_time():
    tr = simulate(**STARDUST_ENTRY, t_end=STARDUST_FLIGHT["drogue_time"], dt_out=0.1)
    assert tr.h[-1] == pytest.approx(STARDUST_FLIGHT["drogue_altitude"], abs=1e3)


def test_conversion_reproduces_pathfinder_published_relative_state():
    from grapes.cases import (MARS_AT_PATHFINDER, PATHFINDER_INERTIAL as s,
                              PATHFINDER_RELATIVE_PUBLISHED as pub)
    v, gam, az = inertial_to_relative(MARS_AT_PATHFINDER, s["h0"], s["v"], s["gamma"],
                                      s["latitude"], s["azimuth"])
    assert v == pytest.approx(pub["v"], abs=1.0)
    assert np.degrees(gam) == pytest.approx(pub["gamma_deg"], abs=0.001)
    assert np.degrees(az) == pytest.approx(pub["azimuth_deg"], abs=0.001)
