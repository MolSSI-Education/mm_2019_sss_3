import numpy as np
import mcpy.particles


class Integrator:
    '''Integrator of a Monte Carlo simulation
    
    Parameters
    ----------
    beta : float
        The inverse of reduced temperature, 1 / T.
    low_acceptance : float, optional,default : 0.38
        The accept rate considered to be low, should be a positive
        float less than 1.
    high_acceptance : float, optional, default : 0.42
        The accept rate considered to be high, should be a positive
        float less than 1, but larger than low_acceptance
    max_displacement : float, optional, default : 0.1
        The initial maximum displacement value in the move. 

    Return
    ------
    self : Integrator
        Return an instance of itself.

    Attributes
    ----------
    beta : float
        The inverse of reduced temperature, 1 / T.
    low_acceptance : float
        The accept rate considered to be low, should be a positive
        float less than 1.
    high_acceptance : float
        The accept rate considered to be high, should be a positive
        float less than 1, but larger than low_acceptance
    max_displacement : float
        The initial maximum displacement value in the move.
    '''
    
    
    def __init__(self,
                 beta,
                 low_acceptance=0.38,
                 high_acceptance=0.42,
                 max_displacement=0.1):
        self.beta = beta
        self.low_acceptance = low_acceptance
        self.high_acceptance = high_acceptance
        self.max_displacement = max_displacement

    def get_particle_energy(self,
                            potential,
                            particles,
                            box_object,
                            i_particle,):
        ''' Computes the energy of a particle with the rest of the system.

        Parameters
        ----------
        potential : class Pairwise object
            A pairwise potential object that can calculate an energy given an
            array of squared distances.
        particles : class object
            particles.coordinates is what will be used. np array.
            An array of atomic coordinates (x, y, z). Shape (n, 3), where
            n is the number of particles.
        box_object: class object
            box_object.box_length is used. float.
            The dimensions of the square box in reduced units.
        i_particle : np.array
            An array of atomic particles (x, y, z). Shape (1, 3).

        Returns
        -------
        e_pair : float
        Total energy of particle i with the rest of the system.
        '''

        rij2 = box_object.minimum_image_distance(i_particle,
                                                 particles.coordinates
                                                 )

        e_pair = potential(rij2)
        # pair_energy_object needs to be replaced by corresponding function.

        return e_pair

    def accept_or_reject(self, delta_e):
        '''Accept or reject a given move based on the Metropolis Criteria.

        Parameters
        ----------
        delta_e : double
            The energy difference between the current step and the previous step.

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

        Parameters
        ----------
        max_displacement : float
            The current value for maximum displacement in the move.
        acc_rate : float
            acceptance rate. It is equal to n_accept/n_trial.

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

    def __call__(self, potential, particles, box, tune_displacement, acc_rate):
        '''Execute a displace trial move.

        Parameters
        ----------
        potential : class Pairwise object
            A pairwise potential object that can calculate an energy given an
            array of squared distances.
        particles : Particles class object
            The Particles class object that holds all infomation of particles.
        box : Box class object
            The Box class object that holds all information of the Box.
        tune_displacement : bool
            If true, integrator tunes the current trial move.
        acc_rate : float
            The acceptance rate for the current trial move
        
        Returns
        -------
        acceptance : bool
            If the current trial move is accepted ( True ), or rejected.
        delta_e : float
            The energy change before and after the current trial move.

        '''
        i_particle = np.random.randint(particles.num_particles)
        random_displacement = (2.0 * np.random.rand(3) - 1.0) * \
            self.max_displacement

        old_energy = self.get_particle_energy(potential,
                                              particles,
                                              box,
                                              i_particle)
        proposed_coordinates = particles.coordinates.copy()
        proposed_coordinates[i_particle] += random_displacement
        proposed_particles = mcpy.particles.Particles(proposed_coordinates)

        new_energy = self.get_particle_energy(potential,
                                              proposed_particles,
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
