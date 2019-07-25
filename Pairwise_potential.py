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
    """Pairwise potential and correction energy by Lennard-Jones potential

    Parameters
    ----------

    sigma : float
        Distance between two particles when interaction potential is zero

    epsilon: float
        Depth of potential well

    """
    
    def __init__(self,sigma=1.0,epsilon=1.0, cutoff=2.6):

        self.sigma = sigma
        self.epsilon = epsilon
        self.cutoff2 = cutoff * cutoff
    
    def lennard_jones(self, rij2):
        """Pairwise potential energy by Lennard-Jones potential

    Parameters
    ----------

    rij2 : float
        Square distance between two particles

    Return
    ------

    Pairwise Lennard-Jones potential energy at considered distance
    
    """
        sig_by_r6 = np.power(self.sigma/rij2,3)
        sig_by_r12 = np.power(sig_by_r6,2)
        return 4.0*self.epsilon*(sig_by_r12-sig_by_r6)

    def cutoff_correction(self, cutoff2, box_object, num_particles,):
        """The function corrects interaction energy from energy cutoff.

    Parameters
    ----------
    cutoff2 : float
        Lennard-Jones potential cutoff distance

    box_object : box
        This is a box object.

    num_particles : float
        Total number of particles in the box

    Return
    ------

    e_correction : float
        Correction energy from truncation

    """

        volume = box_object.volume
        sig_by_cutoff3 = np.power(self.sigma/cutoff2, 3)
        sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
        e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
        e_correction *= 8.0 / 9.0 * np.pi * num_particles / volume * num_particles
        return e_correction

    def __call__(self, rij2):

        if rij2 < self.cutoff2:
            e_pair = lennard_jones(rij2)
        else:
            e_pair = 0.0
        return e_pair

class HS(PairwisePotential):
    """Pairwiswe potential energy by Hard-sphere potential

    Parameters
    ----------

    rij2 = square distance between two particles

    """

    def __call__(self,rij2):
        pass

class SW(PairwisePotential):
    """Pairwiswe potential energy by Square-well potential

    Parameters
    ----------

    rij2 = square distance between two particles

    """

    def __call__(self,rij2):
        pass



