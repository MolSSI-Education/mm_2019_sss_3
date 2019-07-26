# Use call function to create a callables that tries a move accepts or rejects and if accepts,
# update the position, bla bla
import numpy as np


class Integrator:
    def __init__(self,
                 pair_energy_object,
                 cutoff=3.0,
                 low_acceptance=0.38,
                 high_acceptance=0.42):
        self.pair_energy_object = pair_energy_object
        self.cutoff = cutoff
        self.low_acceptance = low_acceptance
        self.high_acceptance = high_acceptance


    def get_particle_energy(self,
                   coordinates,
                   i_particle,
                   box_object):
        ''' Computes the energy of a particle with the rest of the system.

        ----------
        Parameters
        ----------
        coordinates : np.array
        An array of atomic coordinates (x, y, z). Shape (n, 3), where
        n is the number of particles.
        box_length : float
        The dimensions of the square box in reduced units.
        i_particle : np.array
        An array of atomic particles (x, y, z). Shape (1, 3).

        -------
        Returns
        -------
        e_total : float
        Total energy of particle i with the rest of the system.
        '''

        cutoff2 = np.power(self.cutoff, 2)
        e_total = 0.0

        i_position = coordinates[i_particle]

        particle_count = len(coordinates)

        for j_particle in range(particle_count):

            if i_particle != j_particle:

               j_position = coordinates[j_particle]

               rij2 = box_object.minimum_image_distance(i_position, j_position)            ##

               if rij2 < cutoff2:
                   e_pair = self.pair_energy_object(rij2)          ##
                   e_total += e_pair

        return e_total


    def is_accepted(self,
                    delta_e,
                    beta):
        '''Accept or reject a given move based on the Metropolis Criteria.

        ----------
        Parameters
        ----------
        delta_e : double
            The energy difference between the current step and the previous step.
        beta : double
            The inverse of reduced temperature, 1 / T

        -------
        Returns
        -------
        accept : bool
            If the move is accepted (true) or rejected (false).
        '''

        if delta_e < 0.0:
            accept = True

        else:
            random_number = np.random.rand(1)
            p_acc = np.exp(-beta * delta_e)

            if random_number < p_acc:
                accept = True
            else:
                accept = False

        return accept

    def adjust_displacement(self,
                            max_displacement,
                            n_trials,
                            n_accept):
        '''Change max trial displacement on the fly based on acceptance rate.

        ----------
        Parameters
        ----------
        n_trials : int
            The current number of attempted MC moves.
        n_accept : int
            The current number of accepted MC moves.

        -------
        Returns
        -------
        n_trials : int
            Returns zero to reset acceptance rate
        n_accept : int
            Returns zero to reset acceptance rate
        '''

        acc_rate = float(n_accept) / float(n_trials)

        if (acc_rate < self.low_acceptance):
            max_displacement *= 0.8

        elif (acc_rate > self.high_acceptance):
            max_displacement *= 1.2

        n_trials = 0
        n_accept = 0

        return max_displacement, n_trials, n_accept

    @property
    def max_displacement(self):
        max_displacement = .1
        return max_displacement

    @max_displacement.setter
    def new_max_displacement(self, new_max_displacement):
        self.max_displace = self.adjust_displacement(max_displacement, )
        n_trial = 0
        n_accept = 0
        return

    max_displacement = 0.1
    def __call__(self, particles, box):

        i_particle = np.random.randint(particles.num_particles)
        random_displacement = (2.0 * np.random.rand(3) - 1.0) * max_displacement
        cutoff2 = np.power(self.cutoff, 2)   # Are we feeding cutoff here?
        old_energy = self.get_particle_energy(coordinates, box.box_length, i_particle, cutoff2)
        proposed_coordinates = particles.coordinates.copy()
        proposed_coordinates[i_particle] += random_displacement
        acceptance = self.is_accepted()
        new_energy = self.get_particle_energy(new_coordinate)
        delta_e = new_energy - old_energy
        if np.mod(i_step + 1, freq) == 0:

            print(i_step + 1, energy_array[i_step])

            if tune_displacement:
                max_displacement = self.adjust_displacement(n_trials, n_accept, max_displacement)
                n_accept = 0
                n_trials = 0
        return acceptance, delta_e



# try_move = Integrator(coordinates)  ==> this should return accept status. if accepted then return delta e

