"""The planar model against the Allen-Eggers closed-form ballistic entry.

Allen & Eggers, "A Study of the Motion and Aerodynamic Heating of Ballistic
Missiles Entering the Earth's Atmosphere at High Supersonic Speeds",
NACA Report 1381, 1958.

Assumptions of the solution: exponential atmosphere, no gravity, no lift,
flat planet, so the flight-path angle stays at its entry value. Under them:

    V(h)          = V0 exp(-rho(h) H / (2 B sin|gamma|)),  B = m / (Cd A)
    peak decel    = V0^2 sin|gamma| / (2 e H),  where V = V0 exp(-1/2)
    peak heating  (q ~ sqrt(rho) V^3) where V = V0 exp(-1/6)

V0 is the speed far above the atmosphere. The model starts at 150 km, where
a little drag has already acted, so V0 is recovered from the start speed.
The model is run with gravity switched off and a huge planet radius, so it
should reproduce these to integrator precision.
"""
import numpy as np
import pytest

from grapes import Vehicle, simulate
from grapes.atmosphere import Exponential
from grapes.bodies import Body

RHO0, H = 1.225, 7200.0  # sea-level density and scale height for the test
FLAT = Body("flat, gravity-free", radius=1e15, mu=0.0,
            atmosphere=Exponential(RHO0, H), k_sutton_graves=1.75e-4)
VEHICLE = Vehicle("test capsule", mass=46.0, ref_area=0.51887, cd=1.49, nose_radius=0.2286)
B = VEHICLE.ballistic_coefficient


@pytest.fixture(scope="module", params=[-8.2, -30.0, -65.0])
def case(request):
    gamma0 = np.radians(request.param)
    v0, h0 = 12.9e3, 150e3
    t_end = 4 * (h0 - 20e3) / (v0 * np.sin(-gamma0))  # long enough to pass both peaks
    tr = simulate(FLAT, VEHICLE, h0, v0, gamma0, t_end, dt_out=0.002,
                  stop_altitude=10e3, rtol=1e-11, atol=1e-9)
    v_inf = v0 * np.exp(FLAT.density(h0) * H / (2 * B * np.sin(-gamma0)))
    return tr, v_inf, gamma0


def test_flight_path_angle_stays_constant(case):
    tr, _, gamma0 = case
    np.testing.assert_allclose(tr.gamma, gamma0, atol=1e-9)


def test_speed_follows_the_closed_form(case):
    tr, v0, gamma0 = case
    v_exact = v0 * np.exp(-FLAT.density(tr.h) * H / (2 * B * np.sin(-gamma0)))
    np.testing.assert_allclose(tr.v, v_exact, rtol=1e-7)


def test_peak_deceleration(case):
    tr, v0, gamma0 = case
    peak = v0**2 * np.sin(-gamma0) / (2 * np.e * H)
    # 1e-5: the output spacing can miss the exact top of a steep peak
    assert tr.decel.max() == pytest.approx(peak, rel=1e-5)
    assert tr.v[tr.decel.argmax()] == pytest.approx(v0 * np.exp(-0.5), rel=1e-3)


def test_peak_heating_speed(case):
    tr, v0, _ = case
    assert tr.v[tr.q_conv.argmax()] == pytest.approx(v0 * np.exp(-1 / 6), rel=1e-3)
