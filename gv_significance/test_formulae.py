import pytest
import numpy as np

from poisson_poisson import significance as poi_poi_significance
from poisson_poisson import z_bi_significance as poi_poi_zbi

from poisson_gaussian import significance as poi_gau_significance
from poisson_gaussian import z_bi_significance as poi_gau_zbi

from ideal_case import significance as ideal_significance


def test_poi_poi_significance():

    # Test with one number
    s = poi_poi_significance(20, 80, 0.1, k=0, sigma=0)

    s = poi_poi_significance(20, 80, 0.1, k=0.1, sigma=0)

    s = poi_poi_significance(20, 80, 0.1, k=0, sigma=0.1)

    assert np.isclose(poi_poi_zbi(140, 100, 1 / 1.2), 3.93, rtol=0.01)

    # Test with arrays
    s = poi_poi_significance([20, 30], [80, 90], 0.1)
    s = poi_poi_significance([20, 30], [80, 90], [0.1, 0.2])

    s = poi_poi_significance([20, 30], [80, 90], 0.1, k=0.1)
    s = poi_poi_significance([20, 30], [80, 90], [0.1, 0.2], k=[0.1, 0.1])

    s = poi_poi_significance([20, 30], [80, 90], 0.1, sigma=0.1)
    s = poi_poi_significance([20, 30], [80, 90], [0.1, 0.2], sigma=[0.1, 0.1])

    assert np.allclose(poi_poi_zbi([140, 140], [100, 100], [1 / 1.2, 1/1.2]), [3.93, 3.93], rtol=0.01)

    # Test mixed
    s = poi_poi_significance([20, 30, 40], [80, 90, 100], alpha=[0.1, 0.1, 0.1],
                             k=[0, 0.1, 0], sigma=[0, 0, 10])

def test_poi_gau_significance():

    # Test with one number
    s = poi_gau_significance(120, 80, 5.3)

    assert np.isclose(poi_gau_zbi(140, 83.33, 8.333), 3.93, rtol=0.01)

    # Test with vectors
    s = poi_gau_significance([120, 130], [80, 90], [5.3, 5.3])

    assert np.allclose(poi_gau_zbi([140, 140], [83.33, 83.33], [8.333, 8.333]), [3.93, 3.93], rtol=0.01)
