import numpy as np
import pytest

from grapes.atmosphere import Exponential, Table, us_standard_1976

# U.S. Standard Atmosphere 1976, geometric altitude (km) -> density (kg/m^3).
US76 = {0: 1.2250, 11: 3.6480e-1, 20: 8.8910e-2, 32: 1.3555e-2, 50: 1.0269e-3,
        60: 3.0968e-4, 70: 8.2829e-5, 80: 1.8458e-5, 100: 5.604e-7, 120: 2.222e-8}


@pytest.mark.parametrize("km, rho", US76.items())
def test_us76_reference_values(km, rho):
    assert us_standard_1976().density(km * 1e3) == pytest.approx(rho, rel=2e-4)


def test_us76_is_smooth_across_the_86_km_join():
    atm = us_standard_1976()
    below, above = atm.density(85.9e3), atm.density(86.1e3)
    assert above < below < 1.05 * above


def test_us76_always_decreases():
    h = np.linspace(-1e3, 200e3, 2000)
    assert np.all(np.diff(us_standard_1976().density(h)) < 0)


def test_exponential():
    atm = Exponential(rho0=2.0, scale_height=1000.0)
    assert atm.density(1000.0) == pytest.approx(2.0 / np.e)


def test_table_matches_exponential_it_was_built_from():
    exp = Exponential(1.2, 8000.0)
    h = np.linspace(0, 100e3, 11)
    tab = Table(h, exp.density(h))
    for x in (-5e3, 37e3, 140e3):  # below, inside, above the table
        assert tab.density(x) == pytest.approx(exp.density(x))


def test_table_rejects_unsorted_altitudes():
    with pytest.raises(ValueError):
        Table([0, 2, 1], [1, 0.5, 0.7])
