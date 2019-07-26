import numpy as np

class Integrator:
    def __init__(self, 
                pair_energy_object, 
                cutoff = 3.0, 
                low_acceptance = 0.38, 
                high_acceptance = 0.42):
        self.pair_energy_object = pair_energy_object
        self.cutoff = cutoff
        self.low_acceptance = low_acceptance
        self.high_acceptance = high_acceptance

    def __str__(self):
        return 'Intergrator:' 

## tocheck /w box_object & self.pair_energy_object
    def get_particle_energy(self, 
                            coordinates, 
                            i_position, 
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
        i_position : np.array
            An array containing the coordinate of particle i. np.array.shape (1, 3) -> [x, y, z]

        -------
        Returns
        -------
        e_total : float
            Total energy of particle i with the rest of the system.
        '''

        cutoff2 = np.power(self.cutoff, 2)
        e_total = 0.0

        particle_count = len(coordinates)
        
        rij2_array = box_object.minimum_image_distance(i_position, coordinates) 
        
        for rij2 in rij2_array:
        
            if rij2 <= cutoff2:
                e_total += pair_energy_object.potential(rij2)
        
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

