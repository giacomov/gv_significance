import scipy.special
import numpy as np


# This is the smallest number we can process within the numerical precision available
tiny = np.finfo(float).tiny


def significance_from_pvalue(pvalue):
    """
    Return the significance (i.e., the z-score) for a given probability, i.e., the number of standard deviations (sigma)
    corresponding to the probability

    :param pvalue: input probability
    :return: z-score (or significance in units of sigma)
    """

    # Make sure that we can compute this score (i.e., the pvalue is not too small)
    if np.any(pvalue <= tiny):

        raise ArithmeticError("One or more pvalues are too small for a significance computation.")

    # This is equivalent to sqrt(2) * erfinv(1 - 2 * pvalue) but faster

    return -scipy.special.ndtri(pvalue)