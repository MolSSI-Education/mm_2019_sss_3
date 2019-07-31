from mcpy.integrator import Integrator
from mcpy.box import Box
from mcpy.particles import Particles
from mcpy.pairwise import LJ 
import pytest
import sys
import numpy as np

@pytest.mark.parametrize("max_displacement, acc_rate, expected_max_displacement",[
    (0.1, 0.5, 0.11),
    (0.2, 0.3, 0.18),
    (0.3, 0.4, 0.3)
])
def test_adjust_displacement(max_displacement, acc_rate, expected_max_displacement):
    beta = 1 
    inte = Integrator(beta, 0.38, 0.42, max_displacement)
    calculated_max_displacement = inte.adjust_displacement(max_displacement, acc_rate)

    assert np.isclose(calculated_max_displacement, expected_max_displacement)

def test_accept_or_reject():
    delta_e = -np.log(.9)
    trials = 10000
    beta = 1
    inte = Integrator(beta, 0.38, 0.42, 0.1)
    n_acc = 0
    for i in range(trials):
        if inte.accept_or_reject(delta_e) == True:
            n_acc += 1
    p_acc = n_acc / trials 

    assert abs( p_acc - 0.9 ) <= 0.01

