"""
particles.py
A python package for the MolSSI Software Summer School.
Contains Particles class and subclasses RandomParticles and FetchParticles
"""

import numpy as np
from box import *
from abc import ABC, abstractmethod

class Particles(ABC):
    ''' Abstract Particles class
    '''

    def __init__(self):
        '''The constructor of the Particles Class.
        '''
        self._num_particles = 0
        self._coordinates = []

    def __str__(self):
        return( F'Particles Object: {self._num_particles} particles.' )

    @abstractmethod
    def _generate_coordinates(self):
        '''Generates the coordinates of the Particles Object '''
        pass

    @property
    def coordinates(self):
        '''Returns the coordinates of the Particles Object'''
        return self._coordinates

    @property
    def num_particles(self):
        '''Returns the coordinates of the Particles Object'''
        return self._num_particles

    @property
    def box(self):
        '''Returns the Box Object of the Particles Object'''
        return self._box
    
    def get_num_particles(self):
        '''Returns the number of particles of the Particles Object'''
        self._num_particles = len(self._coordinates)


class RandomParticles(Particles):
    ''' This class returns a Particles Object with random positions given number of particles an box dimensions. 
    '''
    def __init__(self, num_particles, box_dims):
        ''' The constructor of the RandomParticles Class.

        Parameters:
        -----------
            num_particles : int
                Number of particles of the system to evaluate. 
            box_dims : numpy.array
                An array of shape (3,) with the dimensions of the box (x, y, z).
        '''
        super().__init__()
        self._num_particles = num_particles
        self._box = Box(box_dims)
        self._coordinates = self._generate_coordinates()
    
    def _generate_coordinates(self):
        '''Generates the coordinates of the RandomParticles Object 

            Returns
            -------
            coordinates : np.array
                Arrays of the wrapped atomic coordinates.
        '''
        coordinates = (0.5 - np.random.rand(self.num_particles, 3)) * self._box.box_dim[np.newaxis,:]
        return coordinates


class FetchParticles(Particles):
    ''' This class returns a Particles Object with random positions given number of particles an box dimensions. 
    '''

    def __init__(self, file_name, wrap_dims = None): #TODO: skiprows=2, usecols=(1, 2, 3)
        ''' FileParticles Class Constructor.

        Parameters:
        -----------
            file_name : str
                A string with the path to the file with x, y, z coordinates.
            wrap_dims : None, np.array
                Determines wheter wraps the particles within the box dimensions (np.array, shape=(3,)) 
                or takes the length dimensions from the file (None). Default None.
        '''
        super().__init__()
        self._file_name = file_name
        self._wrap_dims = wrap_dims
        self._coordinates = self._generate_coordinates()
        self._num_particles = self.get_num_particles()
        self._box = self._create_box()
        
    def _generate_coordinates(self):
        '''Generates the coordinates of the FetchParticles Object 

            Returns
            -------
            coordinates : np.array
                Arrays of the wrapped atomic coordinates.'''
        coordinates = np.loadtxt(self._file_name, skiprows=2, usecols=(1, 2, 3))
        return coordinates

    def get_num_particles(self):
        return (len(self.coordinates))

    def _create_box(self):
        '''Generates the Box object the FetchParticles Object 
        
            Returns
            -------
            box : Box class object
                An Box class object.'''
        if self._wrap_dims is None:
            box_dims = abs( self.coordinates.max(axis=0) - self.coordinates.min(axis=0) )
            box = Box(box_dims)

        elif isinstance(self._wrap_dims, np.ndarray):
            box_dims = self._wrap_dims
            box = Box(box_dims)
            self._coordinates = box.wrap(self._coordinates)

        return box


# Uncoment to test
'''
num_particles = 3
box_dims = np.array([15, 15, 15])
part =  RandomParticles(num_particles, box_dims)
#part = FetchParticles('./sample_config1.xyz', wrap_dims = box_dims)

# Creating a Particles Object
print( part.coordinates)
print(F'\nNum particles: {part.num_particles} \nBox dims: {part.box.box_dim}')
# Rewraping the system example
part.box.box_dim = np.array([1, 1, 1])
part.box.wrap(part.coordinates)
print(part.coordinates)
print(F'\nNum particles: {part.num_particles} \nBox dims: {part.box.box_dim}')
print(part)
'''