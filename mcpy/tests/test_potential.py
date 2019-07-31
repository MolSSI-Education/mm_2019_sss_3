"""
Unit test for the Pairwise_potential calculation.
"""
from mcpy.pairwise import LJ
import pytest
import sys
import numpy as np

def test_pairwise_potential():
    """Three cases are tested:
            - 2 Argon atoms
            - All distances are larger than cutoff
            - Various distances including smaller and larger than cutoff
    """

    # Argon atoms:
    epsilon = 0.997
    sigma = 3.40

    cutoff_ar = 20
    lj_ar_class = LJ(sigma, epsilon, cutoff_ar)

    rij2 = 16
    U_Ar_expected = -0.9368 
    U_Ar_calculated = lj_ar_class.potential(rij2)
    
    # Set of distances that is larger than default cutoff:
    epsilon = 1.0
    sigma = 1.0
    cutoff_large = 2.6
    lj_default_class = LJ(sigma, epsilon, cutoff_large)
    rij2 = np.array([10,10,10,10,10,10,10])
    U_large_expected = 0.0
    U_large_calculated = lj_default_class(rij2)

    # Set of mixed distances:
    rij2 = np.array([10,5,10,5,10,5,10])
    U_mixed_expected = 3*4*(np.power(1/5,6)-np.power(1/5,3))
    U_mixed_calculated = lj_default_class(rij2)

    assert np.isclose(U_Ar_expected,U_Ar_calculated,rtol = 1e-4)
    assert np.isclose(U_large_expected,U_large_calculated)
    assert np.isclose(U_mixed_expected,U_mixed_calculated)


