import numpy as np


def size_one_or_n(value, other_array, name):

    value_ = np.array(value, dtype=float, ndmin=1)

    if value_.shape[0] == 1:

        value_ = np.zeros(other_array.shape[0], dtype=float) + value

    else:

        assert value_.shape[0] == other_array.shape[0], "The size of %s must be either 1 or the same size of n" % name

    return value_
