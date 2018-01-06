import numpy as np
from xlogy import xlogy
import scipy.optimize


def _li_and_ma(n_, b_, alpha):

    # In order to avoid numerical problem with 0 * log(0), we add a tiny number to n_ and b_, inconsequential
    # for the computation
    n_ += 1e-15  # type: np.ndarray
    b_ += 1e-15  # type: np.ndarray

    # Pre-compute for speed
    n_plus_b = n_ + b_
    ap1 = alpha + 1

    res = n_ * np.log(ap1 / alpha * (n_ / n_plus_b))

    res += b_ * np.log(ap1 * (b_ / n_plus_b))

    return np.sqrt(2 * res)


def _likelihood_with_sys(o, b, a, s, k, B, M):
    # Keep away from unphysical solutions during maximization
    # (where the argument of the logarithm is negative) by returning
    # a large negative number

    if M + a * B <= 0 or k + 1 <= 0 or B <= 0:
        return -1000

    # Pre-compute for speed

    Ba = B * a
    Bak = B * a * k

    res = -Bak - Ba - B - M + xlogy(b, B) - k ** 2 / (2 * s ** 2) + xlogy(o, Bak + Ba + M)

    return res


def _get_TS_by_numerical_optimization(n_, b_, alpha, sigma):

    # Numerical optimization to find the maximum likelihood value
    # for the null hypothesis (M=0). We use the closed form for B_mle provided in the paper
    # so that the optimization is in one variable (kk)
    # NOTE: optimize.minimize minimizes, so we multiply the log-likelihood by -1

    wrapper = lambda kk: -1 * _likelihood_with_sys(n_, b_, alpha, sigma, kk,
                                                   B=(b_ + n_) / (alpha * kk + alpha + 1),
                                                   M=0)

    res = scipy.optimize.minimize(wrapper,
                                  [0.0],
                                  tol=1e-3)

    # Get the minimum of the -log likelihood
    h0_mlike_value = res['fun']

    # Get alternative hypothesis minimum of the -log likelihood

    h1_mlike_value = -(xlogy(b_, b_) - b_ + xlogy(n_, n_) - n_)

    # Compute Test Statistic of Likelihood Ratio Test

    TS = 2 * (h0_mlike_value - h1_mlike_value)

    return TS


def significance(n, b, alpha, sigma=0, k=0):
    """
    Returns the significance for detecting n counts when alpha * B are expected.

    If sigma=0 and k=0 (default), this is the case with no additional systematic error and the classic result
    from Li & Ma (1983) is used. Example:

        > significance(n, b, alpha)

    If k>0 then eq.7 from Vianello (2018) is used, which assumes that k is the upper boundary on the fractional
    systematic uncertainty. In this case sigma has no meaning and is ignored. Example:

        > significance(n, b, alpha, k=0.1)

    If sigma>0, then eq. 9 from Vianello (2018) is used, which assumes a Gaussian distribution for the systematic
    uncertainty. In this case k has no meaning and is ignored.Example:

        > significance(n, b, alpha, sigma=0.1)

    :param n: observed counts (can be an array)
    :param b: expected background counts (can be an array)
    :param alpha: ratio of the source observation efficiency and background observation efficiency (must be the same for
    all items in n and b)
    :param sigma: standard deviation for the Gaussian case (must be the same for all items in n and b)
    :param k: maximum fractional systematic uncertainty expected. (must be the same for all items in n and b)
    :return: the significance (z score) for the measurement(s)
    """

    # Make sure we are dealing with arrays, and if not, make the input so. This way
    # we can unify the treatment

    n_ = np.array(n, dtype=float, ndmin=1)
    b_ = np.array(b, dtype=float, ndmin=1)

    # Assign sign depending on whether n_ > b_

    sign = np.where(n_ >= alpha * b_, 1, -1)

    if sigma == 0 and k == 0:

        # Li & Ma
        return np.squeeze(sign * _li_and_ma(n_, b_, alpha))

    if k > 0:

        # Need to use eq. 7 from Vianello (2018), which is simply Li & Ma with alpha -> alpha * (k+1)
        return np.squeeze(sign * _li_and_ma(n_, b_, alpha * (k+1)))

    if sigma > 0:

        # Need to use eq. 9 from Vianello (2018)

        # The computation is not vectorized because there is an optimization involved, so we can only loop...
        TS = np.array([_get_TS_by_numerical_optimization(n_[i], b_[i], alpha, sigma) for i in range(n_.shape[0])])

        # Return significance

        return np.squeeze(sign * np.sqrt(TS))