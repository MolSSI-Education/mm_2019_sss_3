
import numpy as np


class Integrator:
    def __init__(self,
                 pair_energy_object,
                 cutoff=3.0,
                 low_acceptance=0.38,
                 high_acceptance=0.42,
                 max_displacement=0.1):
        self.pair_energy_object = pair_energy_object
        self.cutoff = cutoff
        self.low_acceptance = low_acceptance
        self.high_acceptance = high_acceptance
        self.max_displacement= max_displacement


    def get_particle_energy(self,
                   particles,
                    box_object,
                   i_particle,
                   ):
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

        i_position = particles.coordinates[i_particle]

        rij2 = box_object.minimum_image_distance(i_position,
                      particles.coordinates[i_position != particles.coordinates]
                                                     )
            ##

        e_pair = self.pair_energy_object(rij2)

        return e_pair


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
                            acc_rate):
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


        if (acc_rate < self.low_acceptance):
            max_displacement *= 0.8

        elif (acc_rate > self.high_acceptance):
            max_displacement *= 1.2

        return max_displacement




    def __call__(self, particles, box, beta):

        i_particle = np.random.randint(particles.num_particles)
        random_displacement = (2.0 * np.random.rand(3) - 1.0) * \
                              self.max_displacement

        old_energy = self.get_particle_energy(particles.coordinates,
                                              box.box_length,
                                              i_particle)
        proposed_coordinates = particles.coordinates.copy()
        proposed_coordinates[i_particle] += random_displacement

        new_energy = self.get_particle_energy(new_coordinate)
        delta_e = new_energy - old_energy

        acceptance = self.is_accepted(delta_e, beta)
        if acceptance is True:
            particles.coordinates[i_particle] = proposed_coordinates[i_particle]
        if tune_displacement:
            new_max_displacement = self.adjust_displacement(max_displacement,
                                                            acc_rate)
            self.max_displacement(new_max_displacement)

        return acceptance, delta_e



