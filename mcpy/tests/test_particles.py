from mcpy.particles import Particles
import pytest
import sys
import numpy as np 


def test_from_file():
    try:
        particles =  Particles.from_file(file_name = 'mcpy/tests/sample_config.xyz')
    except IOError:
        pass
    else:
        # File includes 800 particles
        expected_num_particles = 800
        calculated_num_particles = particles.num_particles
        assert (expected_num_particles == calculated_num_particles)


def test_from_random():
    num_particles = 3
    box_dims = np.array([1, 1, 1])
    particles =  Particles.from_random( num_particles = 3, box_dims = box_dims)
    expected_num_particles = 3
    calculated_num_particles = particles.num_particles
    
    # Testing coordinates are inside the box dimensions
    particles_within_the_box = (np.abs(particles.coordinates) <= box_dims).all()

    assert(expected_num_particles == calculated_num_particles)
    assert(particles_within_the_box)
