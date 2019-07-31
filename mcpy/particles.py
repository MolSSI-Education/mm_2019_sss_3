"""
particles.py
A python package for the MolSSI Software Summer School.
Contains Particles class
"""

import numpy as np
from mcpy.box import *

class Particles():
    '''Holds all the information of the particles.

    Parameters
    ----------
    coordinates : np.array
        The numpy.array that holds the coordinates of particles with the 
        shape (n, 3), where n is the number of particles.

    Returns
    -------
    self : Particles
        Returns an instance of itself.

    Attributes
    ----------
    coordinates : np.array
        The numpy.array that holds the coordinates of particles with the
        shape (n, 3), where n is the number of particles.
    '''
    def __init__(self, coordinates):
        ''' Particles Class Constructor.

        Parameters:
        -----------
            coordinates : np.array
                Array of shape (n, 3) where n is the number of particles.
        '''
        self.coordinates = coordinates


    def __str__(self):
        return( F'Particles Object: {self._num_particles} particles.' )

    @classmethod
    def from_file(cls, file_name):
        ''' Class method: generates particles from file.

        Parameters:
        -----------
            file_name : str
                A string with the path to the file with x, y, z coordinates.
        Returns:
        --------
            particles : Particles class object
                Particles class object.
        '''
        coordinates = np.loadtxt(file_name, 
            skiprows=2, usecols=(1, 2, 3))
        particles = cls(coordinates)
        return ( particles )

    @classmethod
    def from_random(cls, num_particles, box_dims):
        ''' Class method: generates particles from file.

        Parameters:
        -----------
            num_particles : int
                Number of particles of the system to generate ramdomly if file_name is not given.
            box_dims : np.array
                Array of shape (3,) of x, y, z dimensions.
        
        Returns:
        --------
            particles : Particles class object
                Particles class object.
        '''
        coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_dims[np.newaxis,:]
        particles = cls(coordinates)
        return ( particles )

    @property
    def num_particles(self):
        '''Returns the coordinates of the Particles Object'''
        return len(self.coordinates)
