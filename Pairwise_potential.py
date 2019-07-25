import numpy as np
from abc import ABC, abstractmethod
class PairwisePotential(ABC):
    """Pairwiswe potential energy at considered distance

    Parameters
    ----------

    rij2 = square distance between two particles

    """
    @abstractmethod
    def __call__(self,rij2):
        pass

class LJ(PairwisePotential):
    
    def __init__(self,sigma=1.0,epsilon=1.0):
        self.sigma = sigma
        self.epsilon = epsilon

    def __call__(self, rij2):
        """Pairwise potential energy by Lennard-Jones potential

    Parameters
    ----------

    sigma : float
        distance between two particles when interaction potential is zero

    epsilon: float
        depth of potential well

    rij2 : float
        square distance between two particles

    Return
    ------

    Value of Lennard-Jones potential energy between two particles at considered distance
    
    """
        sig_by_r6 = np.power(self.sigma/rij2,3)
        sig_by_r12 = np.power(sig_by_r6,2)
        return 4.0*self.epsilon*(sig_by_r12-sig_by_r6)

    def cutoff_correction(self, cutoff, box_object, num_particles,):
        '''The function corrects interaction energy from energy cutoff.

    Parameters
    ----------
    cutoff : float
        Lennard-Jones potential cutoff distance

    box_object : np.array
        the NVT box

    num_particles : float
        total number of particles in the box

    Return
    ------

    e_correction : float
        correction energy from truncation
    '''
        volume = np.power(box_object.length, 3)
        sig_by_cutoff3 = np.power(self.sigma/cutoff, 3)
        sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
        e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
        e_correction *= 8.0 / 9.0 * np.pi * num_particles / volume * num_particles
        return e_correction
class HW(PairwisePotential):

    def __call__(self,rij2):
        pass

class SW(PairwisePotential):
    
    def __call__(self,rij2):
        pass



