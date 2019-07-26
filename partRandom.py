import numpy as np
from box import *
from abc import ABC, abstractmethod

class Particles(ABC):
    ''' Abstract Particles class
    '''

    def __init__(self):
        ''' RandomGenerator Class Constructor.

        Parameters:
        -----------
            num_particles (int) : Number of particles of the system to evaluate. 
            box (numpy.array) : An array of shape (3,) with the dimensions of the box (x, y, z). 
        '''
        self._num_particles = 0
        self._coordinates = []

    def __str__(self):
        return( F'Particles Object: {self._num_particles} particles.' )

    @property
    def coordinates(self):
        return self._coordinates

    @abstractmethod
    def _generate_coordinates(self):
        pass

    @property
    def num_particles(self):
        return self._num_particles
    
    def get_num_particles(self):
        self._num_particles = len(self._coordinates)


class RandomParticles(Particles):
    ''' This class returns a Particles Object with random positions given number of particles an box dimensions. 
    '''
    def __init__(self, num_particles, box_dims):
        ''' RandomParticles Class Constructor.

        Parameters:
        -----------
            num_particles (int) : Number of particles of the system to evaluate. 
            box (numpy.array) : An array of shape (3,) with the dimensions of the box (x, y, z). 
        '''
        super().__init__()
        self._num_particles = num_particles
        self._box = Box(box_dims)
        self._coordinates = self._generate_coordinates(self._num_particles)

    @property
    def coordinates(self):
        return self._coordinates
    
    def _generate_coordinates(self, num_particles):
        coordinates = (0.5 - np.random.rand(num_particles, 3)) * self._box.box_dim[np.newaxis,:]
        return coordinates


class FetchParticles(Particles):
    ''' This class returns a Particles Object with random positions given number of particles an box dimensions. 
    '''

    def __init__(self, file_name, wrap_dims = None): #TODO: skiprows=2, usecols=(1, 2, 3)
        ''' FileParticles Class Constructor.

        Parameters:
        -----------
            file_name (str) : A string with the path to the file with x, y, z coordinates.
        '''
        super().__init__()
        self._file_name = file_name
        self._wrap_dims = wrap_dims
        #self._skiprows = skiprows
        #self._usecols = usecols
        self._coordinates = self._generate_coordinates()
        self._num_particles = self.get_num_particles()
        self._box = self._create_box()
        
    def _generate_coordinates(self):
        coordinates = np.loadtxt(self._file_name, skiprows=2, usecols=(1, 2, 3))
        return coordinates

    def get_num_particles(self):
        return (len(self.coordinates))

    def _create_box(self):
        if self._wrap_dims is None:
            box_dims = abs( self.coordinates.max(axis=0) - self.coordinates.min(axis=0) )
            box = Box(box_dims)

        elif isinstance(self._wrap_dims, np.ndarray):
            box_dims = self._wrap_dims
            box = Box(box_dims).wrap(box_dims)

        return box

num = 30
box_dims = np.array([5, 5, 5])
#part =  RandomParticles(num, box_dims)
part = FetchParticles('./sample_config1.xyz', wrap_dims = box_dims)

print( part.coordinates, part.num_particles, part._box.box_dim)
print(part)