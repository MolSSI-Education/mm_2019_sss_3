import numpy as np


class CounterIndex(object):
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
    tuning : bool
        Flag for whether integrators are tuned for acceptance rate.
    frequency : int
        Determines when to log data to std_out, default 10000.
    '''

    def __init__(self, beta=1, tune_integrators=True, frequency=10000):
        self._tuning = tune_integrators
        self.beta = beta
        self.frequency = frequency
        self.step = 0
        self.steps_accepted = []
        self._integrators = []
        self._log_index = CounterIndex()

    def run(self, steps):
        self.check_state()
        self._initialize_state(steps)
        for i in range(steps):
            for i, integrator in enumerate(self._integrators):
                self.step += 1
                acceptance_rate = self.steps_accepted[i] / self.step
                accepted, delta_e = integrator(self.particles,
                                               self.box,
                                               self.beta,
                                               tune_displacement=self.tune,
                                               acc_rate=acceptance_rate)
                if accepted:
                    self.steps_accepted[i] += 1
                    self.energy += delta_e
                if self.step % self.frequency == 0:
                    self.print_log()
                    self._update_log()
        self.print_log()
        self._update_log()

    def run_upto(self, step):
        if self.step >= self.step:
            return None
        self.run(step - self.step)

    def _update_log(self):
        index = self._log_index()
        self.energies[index] = self.energy
        self.steps[index] = self.step

    def _initialize_state(self, steps):
        log_num = int(steps // self.frequency) + 1
        if self.step == 0:
            log_num += 1
            self.steps = np.zeros(log_num)
            self.energies = np.zeros(log_num)
            self.energy = self.calculate_total_energy()
            self._update_log()
        else:
            self.steps = np.concatenate(
                (self.steps, np.zeros(log_num)), axis=None)
            self.energies = np.concatenate((self.energies, np.zeros(log_num)),
                                           axis=None)

    def add_integrator(self, integrator):
        self._integrators.append(integrator)
        self.steps_accepted.append(0)

    def add_potential(self, potential):
        self.potential = potential

    def add_box(self, box):
        self.box = box

    def add_particles(self, particles):
        self.particles = particles

    def calculate_total_energy(self):
        e_total = 0
        for i in np.arange(self.particles.num_particles):
            rij2 = self.box.minimum_image_distance(
                self.particles.coordinates[i],
                self.particles.coordinates[i+1:]
            )
            e_total += self.potential(rij2)
        return e_total + self.potential.cutoff_correction(
            self.box,
            self.particles.num_particles)

    def check_state(self):
        if self._integrators == list():
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

    def print_log(self):
        format_str = 'Step {}, Energy {}, Acceptance Rates {}'
        accepted_rates = np.array(self.steps_accepted) / self.step
        print(format_str.format(self.step,
                                self.energy / self.particles.num_particles,
                                accepted_rates))
