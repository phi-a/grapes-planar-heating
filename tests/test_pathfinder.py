"""Mars Pathfinder against flight values (Braun et al., AIAA 97-3663)."""
import numpy as np
import pytest

from grapes import simulate
from grapes.atmosphere import mars_pathfinder_1997
from grapes.cases import PATHFINDER, PATHFINDER_ENTRY, PATHFINDER_FLIGHT as F


@pytest.fixture(scope="module")
def traj():
    return simulate(**PATHFINDER_ENTRY, t_end=F["chute_time"], dt_out=0.1)


def test_ballistic_coefficient_matches_the_paper():
    assert PATHFINDER.ballistic_coefficient == pytest.approx(62.4)
    assert PATHFINDER.cd == pytest.approx(1.70, abs=0.01)


def test_peak_deceleration(traj):
    assert traj.peak_decel_g == pytest.approx(F["peak_g"], rel=0.03)
    assert traj.t[traj.decel.argmax()] == pytest.approx(F["peak_g_time"], abs=3.0)


def test_peak_heating(traj):
    assert traj.peak_q_total / 1e4 == pytest.approx(F["peak_q"], rel=0.20)
    assert traj.t[traj.q_total.argmax()] == pytest.approx(F["peak_q_time"], abs=5.0)


def test_radiative_is_minor_at_pathfinder_speed(traj):
    assert traj.peak_q_rad < 0.1 * traj.peak_q_conv


def test_state_at_parachute_deploy(traj):
    assert traj.h[-1] == pytest.approx(F["chute_altitude"], abs=1e3)
    assert traj.v[-1] == pytest.approx(F["chute_speed"], rel=0.10)


def test_atmosphere_table_returns_its_own_rows():
    atm = mars_pathfinder_1997()
    assert atm.density(40.070e3) == pytest.approx(3.4200e-4, rel=1e-6)
    assert atm.density(8.959e3) == pytest.approx(8.0848e-3, rel=1e-6)
    assert atm.density(140.855e3) == pytest.approx(1.6720e-10, rel=1e-6)
    h = np.linspace(0, 160e3, 2000)
    assert np.all(np.diff(atm.density(h)) < 0)
