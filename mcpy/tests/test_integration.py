import pytest
import contextlib
import numpy as np
import mcpy.particles
import mcpy.box
import mcpy.pairwise
import mcpy.integrator
import mcpy.mcsimulation
from timeit import default_timer as timer


@pytest.fixture
def mcsimulation():
    reduced_temperature = 0.9
    reduced_density = 0.9
    num_part = 500
    box_dims = np.full(3, np.cbrt(num_part / reduced_density))
    box = mcpy.box.Box(box_dims=box_dims)
    part = mcpy.particles.Particles.from_random(num_particles=num_part,
                                                box_dims=box.box_dims)
    lj = mcpy.pairwise.LJ(cutoff=3.)
    intg = mcpy.integrator.Integrator(1/reduced_temperature)
    mc = mcpy.mcsimulation.MCSimulation()
    mc.add_integrator(intg)
    mc.add_box(box)
    mc.add_particles(part)
    mc.add_potential(lj)
    return mc


def test_runs(mcsimulation):
    mcsimulation.run(10, supress_output=True)


def test_energy_convergence(mcsimulation):
    t1 = timer()
    mcsimulation.run(1000000, supress_output=True)
    t2 = timer()
    time = t2 - t1
    assert time < 250
    average_energy = np.mean(mcsimulation.energies[
        mcsimulation.steps > 499000]) / mcsimulation.particles.num_particles
    # assert that energy converges to NIST values
    assert np.isclose(average_energy, -6.1773, atol=2e-1)
