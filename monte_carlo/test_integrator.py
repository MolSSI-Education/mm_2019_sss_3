
from .integrator import Integrator
import pytest
import sys
import numpy as np

@pytest.mark.parametrize("max_displacement, n_trials, n_accept, expected_max_displacement",[
    (0.1, 1000, 500, 0.12),
    (0.2, 1000, 300, 0.16),
    (0.3, 1000, 400, 0.3)
])
def test_adjust_displacement(max_displacement, n_trials, n_accept, expected_max_displacement):
    pair_energy_object = 0
    inte = Integrator(pair_energy_object)
    calculated_max_displacement, returned_trials, returned_accept = inte.adjust_displacement(max_displacement, n_trials, n_accept)
    
    assert (returned_accept, returned_trials) == (0, 0)
    assert np.isclose(calculated_max_displacement, expected_max_displacement)

def test_is_accepted():
    delta_e = -np.log(.9)
    pair_energy_object = 0
    inte = Integrator(pair_energy_object)
    beta = 1.0
    n_acc = 0
    for i in range(1001):
        if inte.is_accepted(delta_e, beta) == True:
            n_acc += 1
    p_acc = n_acc / 1000.0

    assert ( p_acc - 0.9 ) <= 0.01

