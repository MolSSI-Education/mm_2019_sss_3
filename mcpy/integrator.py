
import numpy as np
import mcpy.particles


class Integrator:
    def __init__(self,
                 beta,
                 pair_energy_object,
                 low_acceptance=0.38,
                 high_acceptance=0.42,
                 max_displacement=0.1):
        self.pair_energy_object = pair_energy_object
        self.beta = beta
        self.low_acceptance = low_acceptance
        self.high_acceptance = high_acceptance
        self.max_displacement = max_displacement

    def get_particle_energy(self,
                            particles,
                            box_object,
                            i_particle,
                            ):
        ''' Computes the energy of a particle with the rest of the system.

        ----------
        Parameters
        ----------
        particles : class object
            particles.coordinates is what will be used. np array.
            An array of atomic coordinates (x, y, z). Shape (n, 3), where
            n is the number of particles.
        box_object: class object
            box_object.box_length is used. float.
            The dimensions of the square box in reduced units.
        i_particle : np.array
            An array of atomic particles (x, y, z). Shape (1, 3).

        -------
        Returns
        -------
        e_pair : float
        Total energy of particle i with the rest of the system.
        '''

        i_position = particles.coordinates[i_particle]

        rij2 = box_object.minimum_image_distance(i_position,
                                                 particles.coordinates[
                                                     i_position !=
                                                     particles.coordinates]
                                                 )

        e_pair = self.pair_energy_object(rij2)
        # pair_energy_object needs to be replaced by corresponding function.

        return e_pair

    def accept_or_reject(self, delta_e):
        '''Accept or reject a given move based on the Metropolis Criteria.

        ----------
        Parameters
        ----------
        delta_e : double
            The energy difference between the current step and the previous step.

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
            p_acc = np.exp(-self.beta * delta_e)

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
        max_displacement : float
            The current value for maximum displacement in the move.
        acc_rate : float
            acceptance rate. It is equal to n_accept/n_trial.

        -------
        Returns
        -------
        max_displacement :float
            New max_displacement adjusted based on acc_rate.
        '''

        if acc_rate < self.low_acceptance:
            max_displacement *= 0.9

        elif acc_rate > self.high_acceptance:
            max_displacement *= 1.1

        return max_displacement

    def __call__(self, particles, box, tune_displacement, acc_rate):

        i_particle = np.random.randint(particles.num_particles)
        random_displacement = (2.0 * np.random.rand(3) - 1.0) * \
            self.max_displacement

        old_energy = self.get_particle_energy(particles,
                                              box,
                                              i_particle)
        proposed_coordinates = particles.coordinates.copy()
        proposed_coordinates[i_particle] += random_displacement
        proposed_particles = mcpy.particles.Particles(proposed_coordinates)

        new_energy = self.get_particle_energy(proposed_particles,
                                              box,
                                              i_particle)
        delta_e = new_energy - old_energy

        acceptance = self.accept_or_reject(delta_e)
        if acceptance is True:
            particles.coordinates[i_particle] = \
                proposed_coordinates[i_particle]
        if tune_displacement:
            self.max_displacement = self.adjust_displacement(
                self.max_displacement,
                acc_rate)

        return acceptance, delta_e
