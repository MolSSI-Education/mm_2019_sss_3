from ..integrator import Integrator
from ..box import Box
from ..particles import Particles
from ..pairwise import LJ 
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
    pair_energy_object = 0
    inte = Integrator(beta, pair_energy_object)
    calculated_max_displacement = inte.adjust_displacement(max_displacement, acc_rate)

    assert np.isclose(calculated_max_displacement, expected_max_displacement)

def test_accept_or_reject():
    delta_e = -np.log(.9)
    pair_energy_object = 0
    beta = 1
    inte = Integrator(beta, pair_energy_object) 
    n_acc = 0
    for i in range(1001):
        if inte.accept_or_reject(delta_e) == True:
            n_acc += 1
    p_acc = n_acc / 1000.0

    assert abs( p_acc - 0.9 ) <= 0.01

#def test_get_particle_energy():
#    particles_object = Particles([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 7.0]])
#    box_object = Box([10.0, 10.0, 10.03])
#    beta = 1
#    pair_energy_object = LJ()
#    inte = Integrator(beta,  pair_energy_object)
#    calculated_particle_energy = inte.get_particle_energy(particles_object, box_object, 0) 
#    expected_particle_energy = pair_energy_object.potential(1.0) + pair_energy_object.potential(5.0) + pair_energy_object.potential(7.0)
#    assert np.isclose(calculated_particle_energy, expected_particle_energy)
