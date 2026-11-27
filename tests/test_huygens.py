"""Huygens at Titan against flight values and the pre-flight heating study.

Flight: Striepe et al., AAS 07-226, and the HASI archive.
Heating: Johnston, Hollis & Sutton, AIAA 2006-3426.
"""
from dataclasses import replace

import numpy as np
import pytest

from grapes import TITAN, simulate
from grapes.atmosphere import titan_huygens_2005
from grapes.cases import HUYGENS_ENTRY, HUYGENS_FLIGHT as F, HUYGENS_INERTIAL
from grapes.heating import (HUYGENS_DESIGN_POINTS as P, TITAN_HUYGENS_BOLTZMANN,
                            TITAN_HUYGENS_NONBOLTZMANN, sutton_graves)

RHO, V = np.array(P["rho"]), np.array(P["v"])


@pytest.fixture(scope="module")
def traj():
    return simulate(**HUYGENS_ENTRY, t_end=F["late_time"], dt_out=0.1)


def test_entry_state_from_the_state_vector():
    assert HUYGENS_INERTIAL["h0"] == pytest.approx(1270e3, abs=100.0)
    assert HUYGENS_INERTIAL["v"] == pytest.approx(6031.0, abs=1.0)
    assert np.degrees(HUYGENS_INERTIAL["gamma"]) == pytest.approx(-65.5, abs=0.1)


def test_peak_deceleration(traj):
    assert traj.decel.max() == pytest.approx(F["peak_decel"], rel=0.04)
    assert traj.t[traj.decel.argmax()] == pytest.approx(F["peak_decel_time"], abs=2.0)


def test_altitude_at_last_entry_record(traj):
    assert traj.h[-1] == pytest.approx(F["late_altitude"], abs=4e3)


def test_convective_near_the_design_study():
    q = sutton_graves(RHO, V, 1.25, TITAN) / 1e4
    np.testing.assert_allclose(q, P["q_conv"], rtol=0.10)


@pytest.mark.parametrize("fit, key", [(TITAN_HUYGENS_NONBOLTZMANN, "q_rad_nonboltzmann"),
                                      (TITAN_HUYGENS_BOLTZMANN, "q_rad_boltzmann")])
def test_radiative_fit_reproduces_its_points(fit, key):
    np.testing.assert_allclose(fit.rate(RHO, V) / 1e4, P[key], rtol=1e-9)


def test_radiative_fit_gives_no_data_outside_its_speed_range():
    fit = TITAN_HUYGENS_NONBOLTZMANN
    assert np.isnan(fit.rate(3e-4, 4000.0))
    assert np.isnan(fit.rate(1e-4, 6000.0))
    assert np.isfinite(fit.rate(3e-4, 5000.0))


def test_radiative_rivals_convective_at_titan(traj):
    upper = simulate(**{**HUYGENS_ENTRY, "body": replace(TITAN, radiative=TITAN_HUYGENS_BOLTZMANN)},
                     t_end=F["late_time"], dt_out=0.1)
    assert traj.peak_q_rad > 0.5 * traj.peak_q_conv
    assert upper.peak_q_rad > traj.peak_q_conv


def test_totals_ignore_missing_radiative_data(traj):
    assert np.isnan(traj.q_rad).any()
    assert np.isfinite(traj.q_total).all()
    assert np.isfinite(traj.heat_load_rad).all()


def test_atmosphere_table():
    atm = titan_huygens_2005()
    assert atm.density(0.003e3) == pytest.approx(5.3446, rel=1e-6)
    assert atm.density(100.003e3) == pytest.approx(0.022699, rel=1e-6)
    # above 157 km density comes from pressure and temperature
    assert atm.density(269.7708e3) == pytest.approx(18.2881 * 28.0 / (8314.46 * 186.361), rel=1e-6)
    h = np.linspace(0, 1400e3, 5000)
    assert np.all(np.diff(atm.density(h)) < 0)
