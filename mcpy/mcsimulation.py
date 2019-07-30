'''mcsimulation.py - Holds the MCSimulation class which interfaces with the
other classes in mcpy to provide a simple interface to run MC simulations.
'''

import numpy as np


class CounterIndex(object):
    '''Just a simple counter class that increments when called.'''

    def __init__(self):
        self._cnt = 0

    def __call__(self):
        cnt = self._cnt
        self._cnt += 1
        return cnt


class MCSimulation(object):
    '''Holds all the information to run a Monte Carlo simulation.

    Parameters
    ----------
    tune_integrators : bool, optional
        If true asks integrators to tune the trial move.
    frequency : int, optional
        Determines when to log data to std_out, default 10000.

    Returns
    -------
    self : MCSimulation
        Returns an instance of itself.

    Attributes
    ----------
    particles : Particles
        The particles in the MC simulation. Particles can be wrapped or
        unwrapped.
    box : Box
        The the periodic box the simulation is running in.
    integrators : list of Integrator objects
        A list of  Integrator objects that generate and accept or reject trial
        moves.
    potential : a pairwise potential object
        A list of potential objects that calculate pair potentials given a
        squared distance.
    tune : bool
        Flag for whether integrators are tuned for acceptance rate. Implemented
        so that it returns false on steps that
        `self.step % self.frequency != 0` regardles of initialization.
    frequency : int
        Determines when to log data to std_out, default 10000.
    '''

    def __init__(self, tune_integrators=True, frequency=10000):
        self._tuning = tune_integrators
        self.frequency = frequency
        self.step = 0
        self.steps_accepted = []
        self.integrators = []
        self._log_index = CounterIndex()

    def run(self, steps, supress_output=False):
        '''Runs the simulation for `steps` steps.

        Also logs data for each `self.frequency` steps. Prints log data to
        standard output at `self.frequency` intervals as well. The last step is
        always printed as well.

        Parameters
        ----------
        steps : int
            The number of MC steps to run.

        Return
        ------
        None
        '''
        self.check_state()
        self._initialize_state(steps)
        for i in range(steps):
            for i, integrator in enumerate(self.integrators):
                self.step += 1
                acceptance_rate = self.steps_accepted[i] / self.step
                accepted, delta_e = integrator(self.potential,
                                               self.particles,
                                               self.box,
                                               acc_rate=acceptance_rate,
                                               tune_displacement=self.tune)
                if accepted:
                    self.steps_accepted[i] += 1
                    self.energy += delta_e
                if self.step % self.frequency == 0:
                    self.print_log(supress_output)
                    self._update_log()

    def run_upto(self, step):
        '''Run until the simulation has reached the specified step.

        Does nothing if step has already been exceeded.

        Parameters
        ----------
        step : int
            The MC step to end on.

        Returns
        -------
        None
        '''
        if self.step >= step:
            return None
        self.run(step - self.step)

    def _update_log(self):
        index = self._log_index()
        self.energies[index] = self.energy
        self.steps[index] = self.step

    def _initialize_state(self, steps):
        log_num = steps // self.frequency + 1
        if self.step == 0:
            self.steps = np.zeros(log_num)
            self.energies = np.zeros(log_num)
            self.energy = self.calculate_total_energy()
            self._update_log()
        else:
            self.steps = np.concatenate((self.steps, np.zeros(log_num)),
                                        axis=None)
            self.energies = np.concatenate((self.energies, np.zeros(log_num)),
                                           axis=None)

    def add_integrator(self, integrator):
        '''Add integrator to list of simulation integrators.

        Parameters
        ----------
        integrator : mcpy.Integrator object
            A initialized mcpy.Integrator object that can generate and realize
            MC trial moves.

        Returns
        -------
        None
        '''
        self.integrators.append(integrator)
        self.steps_accepted.append(0)

    def add_potential(self, potential):
        '''Add a pairwise potential to the simulation.

        Parameters
        ----------
        potential : mcpy.Pairwise.Potential object
            A pairwise potential to use in calculating the energy.

        Returns
        -------
        None
        '''
        self.potential = potential

    def add_box(self, box):
        '''Add a simulation box to the simulation.

        Parameters
        ----------
        potential : mcpy.Box object
            A simulation box for the simulation.

        Returns
        -------
        None
        '''
        self.box = box

    def add_particles(self, particles):
        '''Add particles to the simulation.

        Parameters
        ----------
        potential : mcpy.Particles object
            Particles to use for the simulation.

        Returns
        -------
        None
        '''
        self.particles = particles

    def calculate_total_energy(self):
        '''Calculate the current total energy.

        Uses the potential, particles, and box objects. Usually only needs to
        be done at initialization though can be called at any time.
        '''
        e_total = 0
        for i in np.arange(self.particles.num_particles - 1):
            rij2 = self.box.minimum_image_distance(0,
                    self.particles.coordinates[i:]
            )
            e_total += self.potential(rij2)
        return e_total + self.potential.cutoff_correction(
            self.box,
            self.particles.num_particles)

    def check_state(self):
        '''Raises a RuntimeError if self is not ready to run.'''
        if self.integrators == list():
            raise RuntimeError("No integrator defined.")
        if not hasattr(self, 'potential'):
            raise RuntimeError("No potential defined.")
        if not hasattr(self, 'particles'):
            raise RuntimeError("No particles defined.")
        if not hasattr(self, 'box'):
            raise RuntimeError("No box defined.")

    @property
    def tune(self):
        return self.step % self.frequency == 0 if self._tuning else False

    def print_log(self, supress_output=False):
        '''Print out current step, energy, and integrator acceptance rates.'''
        format_str = 'Step {}, Energy {}, Acceptance Rates {}'
        accepted_rates = np.array(self.steps_accepted) / self.step
        if not supress_output:
            print(format_str.format(self.step,
                                    self.energy / self.particles.num_particles,
                                    accepted_rates))
