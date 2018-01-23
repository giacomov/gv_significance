import numpy as np
from numpy import sqrt, squeeze
from xlogy import xlogyv
from z_bi import z_bi_vectorized
from size_one_or_n import size_one_or_n


def significance(n, b, sigma):
    """
    Returns the significance for observing n counts when b are expected. The measurement "b +/- sigma" is returned
    by some kind of background estimation procedure, and b is assumed to be a Gaussian random variable.

    :param n: observed counts
    :param b: estimation of the background coming from some method
    :param sigma: error on the estimation of the background
    :return: the significance of the measurement(s)
    """

    n_ = np.array(n, dtype=float, ndmin=1)
    b_ = np.array(b, dtype=float, ndmin=1)

    sigma_ = size_one_or_n(sigma, n_, "sigma")

    # Assign sign depending on whether n_ > b_

    sign = np.where(n_ >= b_, 1, -1)

    B0_mle = 0.5 * (b_ - sigma_ ** 2 + sqrt(b_ ** 2 - 2 * b_ * sigma_ ** 2 + 4 * n_ * sigma_ ** 2 + sigma_ ** 4))

    # B0_mle could be slightly < 0 (even though it shouldn't) because of the
    # limited numerical precision of the calculator. let's accept as negative as 0.01, and clip
    # at zero to avoid giving results difficult to interpret
    assert np.all(B0_mle > -0.01), "This is a bug. B0_mle cannot be negative."

    B0_mle = np.clip(B0_mle, 0, None)

    return squeeze(sqrt(2) * sqrt(xlogyv(n_, n_ / B0_mle) + (b_ - B0_mle)**2 / (2 * sigma_**2) + B0_mle - n_) * sign)


def z_bi_significance(n, b, sigma):
    """
    Use the estimator Z_Bi from Cousins et al. 2008 to compute the significance

    :param n: observed counts (can be an array)
    :param b: expected background counts (can be an array)
    :param sigma: error on the measurement of b
    :return: the significance (z score) for the measurement(s)
    """

    n_ = np.array(n, dtype=float, ndmin=1)
    b_ = np.array(b, dtype=float, ndmin=1)
    sigma_ = size_one_or_n(sigma, n_, "sigma")

    res = np.zeros(n_.shape[0], dtype=float)

    idx = (sigma_ != 0)

    # When sigma is zero there is no computation possible
    res[~idx] = np.nan

    tau = b_[idx] / sigma_[idx]**2  # Eq. 8 of Cousins et al. 2008
    n_off = b_[idx] * tau

    res[idx] = z_bi_vectorized(n_[idx], n_off, 1.0 / np.abs(tau))

    return np.squeeze(res)
