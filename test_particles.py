from .particles import Particles
import pytest
import sys
import numpy as np 


def test_from_file():
    particles =  Particles.from_file(file_name='sample_config1.xyz')
    # File includes 800 particles
    expected_num_particles = 800
    calculated_num_particles = particles.num_particles
    assert (expected_num_particles == calculated_num_particles)

def test_from_random():
    num_particles = 3
    box_dims = np.array([1, 1, 1])
    particles =  Particles.from_random( num_particles = 3, box_dims = box_dims)
    # File includes 800 particles
    expected_num_particles = 3
    calculated_num_particles = particles.num_particles
    assert (expected_num_particles == calculated_num_particles)
