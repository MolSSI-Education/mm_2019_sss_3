import pytest
import numpy as np
import mcpy.particles
import mcpy.box
import mcpy.pairwise
import mcpy.integrator
import mcpy.mcsimulation
from timeit import default_timer as timer

def test_energy_convergence():
    reduced_temperature = 0.9
    reduced_density = 0.9
    num_part = 500
    box_dims = np.full(3, np.cbrt(num_part / reduced_density))
    box = mcpy.box.Box(box_dims=box_dims)
    part = mcpy.particles.Particles.from_random(num_particles=num_part,
                                                box_dims=box.box_dims)
    lj = mcpy.pairwise.LJ(cutoff=3.)
    intg = mcpy.integrator.Integrator(1/reduced_temperature, lj)
    mc = mcpy.mcsimulation.MCSimulation()
    mc.add_integrator(intg)
    mc.add_box(box)
    mc.add_particles(part)
    mc.add_potential(lj)
    t1 = timer()
    mc.run(1000000)
    t2 = timer()
    average_energy = np.mean(mc.energies[mc.steps > 499000])
    time = t2 - t1
    # assert 1 million steps occur within 250 seconds
    assert time < 250
    # assert that energy converges to NIST values
    assert np.isclose(average_energy, 6.1773, atol=1e-2)
