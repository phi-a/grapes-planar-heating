import numpy as np
import pytest

from grapes import EARTH, MARS, simulate
from grapes.cases import STARDUST_ENTRY
from grapes.heating import TAUBER_SUTTON_EARTH, TAUBER_SUTTON_MARS


def test_earth_hand_calculation():
    # Rn = 1 m so the exponent drops out: 4.736e4 * (1e-4)^1.22 * f(11 km/s = 151)
    q = TAUBER_SUTTON_EARTH.rate(1e-4, 11000.0, 1.0) / 1e4
    assert q == pytest.approx(4.736e4 * 10**-4.88 * 151, rel=1e-9)
    assert q == pytest.approx(94.3, abs=0.1)


def test_mars_hand_calculation():
    q = TAUBER_SUTTON_MARS.rate(1e-4, 7000.0, 1.0) / 1e4
    assert q == pytest.approx(2.35e4 * 10**-4.76 * 8.1, rel=1e-9)
    assert q == pytest.approx(3.31, abs=0.01)


@pytest.mark.parametrize("model", [TAUBER_SUTTON_EARTH, TAUBER_SUTTON_MARS])
def test_zero_below_the_fit_and_rising_inside_it(model):
    lo, hi = model.v_range
    assert model.rate(1e-4, lo - 1.0, 0.5) == 0.0
    q = model.rate(1e-4, np.linspace(lo, hi, 50), 0.5)
    assert np.all(np.diff(q) > 0)


def test_bigger_nose_radiates_more_but_convects_less():
    small = simulate(**STARDUST_ENTRY, t_end=140)
    big_vehicle = type(small.vehicle)("big", 46.0, small.vehicle.ref_area, 1.368, 0.5)
    big = simulate(**{**STARDUST_ENTRY, "vehicle": big_vehicle}, t_end=140)
    assert big.peak_q_rad > small.peak_q_rad
    assert big.peak_q_conv < small.peak_q_conv


def test_stardust_radiative_is_a_modest_share():
    tr = simulate(**STARDUST_ENTRY, t_end=140, dt_out=0.1)
    share = tr.peak_q_rad / tr.peak_q_total
    assert 0.05 < share < 0.30


def test_sutton_graves_mixing_rule():
    from grapes.heating import sutton_graves_coefficient
    # Titan, as used in the code
    assert sutton_graves_coefficient(dict(N2=0.9699, CH4=0.0230, Ar=0.0071)) == pytest.approx(
        1.7403e-4, rel=1e-3)
    # Mars 97% CO2 / 3% N2 reproduces the value already used for Mars
    assert sutton_graves_coefficient(dict(CO2=0.97, N2=0.03)) == pytest.approx(MARS.k_sutton_graves,
                                                                               rel=2e-3)
    # a pure gas gives its own Table II value
    assert sutton_graves_coefficient(dict(N2=1.0)) == pytest.approx(
        0.1112 / (2 * 101325 ** 0.5), rel=1e-12)


def test_bodies_carry_their_model():
    assert EARTH.radiative is TAUBER_SUTTON_EARTH
    assert MARS.radiative is TAUBER_SUTTON_MARS


def test_table_is_interpolated_linearly():
    # Halfway between 11,000 (151) and 11,500 (238) m/s, Rn = 1 m.
    q = TAUBER_SUTTON_EARTH.rate(1e-4, 11250.0, 1.0) / 1e4
    assert q == pytest.approx(4.736e4 * 10**-4.88 * 0.5 * (151 + 238), rel=1e-9)


def test_earth_exponent_caps_follow_the_paper():
    from grapes.heating import _earth_exponent
    rho, v = 1e-5, 10000.0  # low density, low speed: uncapped a is large
    assert _earth_exponent(rho, v, 0.5) == 1.0
    assert _earth_exponent(rho, v, 1.0) == 0.6
    assert _earth_exponent(rho, v, 2.0) == 0.6
    assert _earth_exponent(rho, v, 2.5) == 0.5


def test_radiative_share_grows_with_nose_radius():
    from dataclasses import replace
    from grapes.cases import PATHFINDER_ENTRY

    def share(entry, rn, t_end):
        tr = simulate(**{**entry, "vehicle": replace(entry["vehicle"], nose_radius=rn)},
                      t_end=t_end)
        return tr.peak_q_rad / (tr.peak_q_rad + tr.peak_q_conv)

    earth = [share(STARDUST_ENTRY, r, 140) for r in (0.3, 1.0, 3.0)]
    mars = [share(PATHFINDER_ENTRY, r, 171) for r in (1.0, 5.0, 23.0)]
    assert earth == sorted(earth) and mars == sorted(mars)
    assert earth[0] < 0.5 < earth[2]  # crossover inside the Earth fit's range
    assert mars[0] < 0.5 < mars[2]


def test_at_huygens_size_titan_matches_earth_at_half_the_speed():
    from dataclasses import replace
    from grapes import TITAN
    from grapes.cases import HUYGENS_ENTRY, PATHFINDER_ENTRY
    from grapes.heating import TITAN_HUYGENS_BOLTZMANN

    def share(entry, t_end, body=None):
        run = {**entry, "vehicle": replace(entry["vehicle"], nose_radius=1.25)}
        if body is not None:
            run["body"] = body
        tr = simulate(**run, t_end=t_end)
        return tr.peak_q_rad / (tr.peak_q_rad + tr.peak_q_conv)

    earth = share(STARDUST_ENTRY, 140)
    mars = share(PATHFINDER_ENTRY, 171)
    titan_low = share(HUYGENS_ENTRY, 267.7)
    titan_high = share(HUYGENS_ENTRY, 267.7, replace(TITAN, radiative=TITAN_HUYGENS_BOLTZMANN))
    assert earth == pytest.approx(0.47, abs=0.01)
    assert mars == pytest.approx(0.08, abs=0.01)
    assert titan_low >= earth - 0.005 and titan_high > earth
