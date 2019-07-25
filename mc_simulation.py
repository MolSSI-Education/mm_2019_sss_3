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

    def __init__(self, tune_integrators=True, frequency=10000):
        self.tuning = tune_integrators
        self.frequency = frequency
        self.step = 0
        self.steps_accepted = []
        self._log_index = CounterIndex()

    def run(self, steps):
        self._initialize_state(steps)
        for i in range(steps):
            for i, integrator in enumerate(self.integrators):
                self.steps += 1
                acceptance_rate = self.steps_accepted[i] / self.step
                delta_e, accepted = integrator(self.particles,
                                               self.box,
                                               acceptance_rate,
                                               tune=self.tune())
                if accepted:
                    self.steps_accepted[i] += 1
                self.energy += delta_e
                if self.steps % self.frequency == 0:
                    self.print_log()
                    index = self._log_index
                    self.energies[index] = self.energy
                    self.steps[index] = self.step

    def _initialize_state(self, steps):
        log_num = steps // self.frequency + 1
        if self.step == 0:
            self.steps = np.zeros(log_num)
            self.energies = np.zeros(log_num)
            self.energy = self.calculate_total_energy()
            index = self._log_index()
            self.energies[index] = self.energy
            self.steps[index] = 0
        else:
            self.steps = np.concatenate(self.steps, np.zeros(log_num))
            self.energies = np.concatenate(self.energies, np.zeros(log_num))

    def add_integrator(self, integrator):
        self.integrators.append(integrator)
        self.steps_accepted.append(0)

    def tune(self):
        return self.steps % self.frequency == 0 if self.tuning else False

    def add_box(self, box):
        self.box = box

    def add_particles(self, particles):
        self.particles = particles

    def print_log(self):
        format_str = 'Step {}, Energy {}, Acceptance Rates {}'
        accepted_rates = np.array(self.steps_accepted) / self.steps
        print(format_str.format(self.step,
                                self.energy / self.particles,
                                accepted_rates))
