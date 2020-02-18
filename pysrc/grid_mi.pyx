import cython
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector

cdef extern from "range_mi/grid_mi.hpp"  namespace "range_mi":

    cppclass GridMI:

        GridMI(unsigned int height, unsigned int width)

        void compute_mi(const double * const vacancy, unsigned int num_beams)

        const vector[double] & mi()

"""
This function turns a map of occupancy
probabilities into a map indicating where
the highest mutual information is.
"""
def grid_mi(
        np.ndarray[double, ndim=2, mode="c"] vacancy,
        unsigned int num_beams):

    # Create a mutual information computation device
    mi_computer = new GridMI(vacancy.shape[0], vacancy.shape[1]);

    # Compute the mutual information
    mi_computer.compute_mi(<double *> vacancy.data, num_beams)

    # Copy the output
    return np.array(mi_computer.mi()).reshape((vacancy.shape[0], vacancy.shape[1]))
