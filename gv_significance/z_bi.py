import numpy as np
from ncephes.cprob import incbet
from numba import vectorize, float64
from significance_from_pvalue import significance_from_pvalue


# This decorator vectorize the function for fast execution
@vectorize([float64(float64, float64, float64)])
def z_bi_cephes(n_on, n_off, alpha):

    tau = 1.0 / alpha

    aa = n_on
    bb = n_off + 1
    xx = 1.0 / (1+tau)

    # Checks to avoid Nan in some cases
    if aa <= 0.0 or bb <= 0.0:

        return 0.0

    if xx <= 0.0:

        return 0.0

    if xx >= 1.0:

        return 1.0

    # I use the incbet from cephes instead of the scipy.special.betainc function because the latter has numerical
    # problems in some instances and return Nans, while the incbet from Cephes is more robust

    P_Bi = incbet(aa, bb, xx)

    return significance_from_pvalue(P_Bi)


def z_bi_vectorized(n, b, alpha):
    """
    Use the estimator Z_Bi from Cousins et al. 2008 to compute the significance

    :param n: observed counts (can be an array)
    :param b: expected background counts (can be an array)
    :param alpha: ratio of the source observation efficiency and background observation efficiency (must be the same for
    all items in n and b)
    :return: the significance (z score) for the measurement(s)
    """

    n_ = np.array(n, dtype=float, ndmin=1)
    b_ = np.array(b, dtype=float, ndmin=1)

    assert n_.shape[0] == b_.shape[0], "n and b must have the same size"

    alpha_ = np.array(alpha, dtype=float, ndmin=1)

    if alpha_.shape[0] == 1:

        alpha_ = np.array([alpha] * n_.shape[0])

    else:

        assert alpha_.shape[0] == n_.shape[0], "Alpha must be either a scalar or an array of the same length of n"

    # Assign sign depending on whether n_ > b_

    sign = np.where(n_ >= alpha * b_, 1, -1)

    res = z_bi_cephes(n_, b_, alpha_)

    return np.squeeze(sign * res)