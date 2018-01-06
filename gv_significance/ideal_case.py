import scipy.special
import numpy as np
from significance_from_pvalue import significance_from_pvalue


def significance(n, B):
    """
    Significance of the measurement of n counts when B are expected, for the ideal case where there is no uncertainty
    on the background B.

    See section 3.1 in Vianello (2018)

    :param n: observed counts
    :param B: expected background (assumed with no uncertainty)
    :return: significance (z-score) for the measurement(s)
    """

    # pdtrc is the summation of the Poisson distribution with average B from n to infinity

    pvalue = scipy.special.pdtrc(n, B)

    # Convert to significance (z score) and return

    return significance_from_pvalue(pvalue)


def five_sigma_threshold(B, detection_efficiency):
    """
    Returns the counts that a source must generate in order to be detected at 5 sigma given a background
    (assumed with no uncertainty), for a given efficiency of detection (1 - type II error).

    NOTE: the returned counts are what the source should generate (on top of the background)

    See eq. 3 in Vianello 2018

    :param B: expected background (no uncertainty)
    :param detection_efficiency: one of 50, 90 or 99. Represents the desired detection efficiency, i.e., 1 - p_II where
    p_II is the probability of a type II error
    :return: counts that a source must generate in order to be detected above 5 sigma with the given efficiency
    """

    legal_efficiencies = [50, 90, 99]

    assert detection_efficiency in legal_efficiencies, \
        "Detection efficiency must be one of %s" % ",".join(legal_efficiencies)

    a = [4.053, 7.391, 11.090]
    b = [5.038, 6.356, 7.415]

    idx = legal_efficiencies.index(detection_efficiency)

    return a[idx] + b[idx] * np.sqrt(B)