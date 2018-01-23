from math import log
import numpy as np


def xlogy(x, y):
    """
    This function implements x * log(y) so that if both x and y are 0, the results is zero (and not infinite or nan
    as the computer would return otherwise).

    NOTE: x and y must be numbers, not arrays
    """

    if x == 0.0:

        return 0.0

    else:

        return x * log(y)


def xlogyv(x, y):
    """
    This function implements x * log(y) so that if both x and y are 0, the results is zero (and not infinite or nan
    as the computer would return otherwise).

    Version which accepts numpy.array as inputs.
    """

    x = np.array(x, ndmin=1)
    y = np.array(y, ndmin=1)

    results = np.zeros_like(y)

    idx = (x != 0)

    results[idx] = x[idx] * np.log(y[idx])

    return np.squeeze(results)
