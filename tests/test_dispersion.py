import numpy as np

from grapes import simulate
from grapes.cases import PATHFINDER_ENTRY, PATHFINDER_FLIGHT
from grapes.dispersion import SIGMA, disperse

T = PATHFINDER_FLIGHT["chute_time"]


def test_first_run_is_the_nominal():
    peak, alt = disperse(PATHFINDER_ENTRY, T, n=3)
    nominal = simulate(**PATHFINDER_ENTRY, t_end=T, dt_out=0.1, rtol=1e-8)
    assert peak[0] == nominal.decel.max()
    assert alt[0] == nominal.h[-1]


def test_zero_uncertainty_gives_no_spread():
    peak, alt = disperse(PATHFINDER_ENTRY, T, n=4, sigma={k: 0.0 for k in SIGMA})
    assert np.ptp(peak) == 0.0 and np.ptp(alt) == 0.0


def test_same_seed_same_result():
    a = disperse(PATHFINDER_ENTRY, T, n=6, seed=1)
    b = disperse(PATHFINDER_ENTRY, T, n=6, seed=1)
    c = disperse(PATHFINDER_ENTRY, T, n=6, seed=2)
    np.testing.assert_array_equal(a[0], b[0])
    assert not np.array_equal(a[0], c[0])
