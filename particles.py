import numpy as np
from abc import ABC, abstractmethod

class ParticlesGenerator(ABC):
    ''' Interface ParticlesGenerator to obtain a Particles Object
    '''
    def __init__(self):
        ''' RandomGenerator Class Constructor.

        Parameters:
        -----------
            num_particles (int) : Number of particles of the system to evaluate. 
            box (numpy.array) : An array of shape (3,) with the dimensions of the box (x, y, z). 
        '''
        self._num_particles = num_particles
        self._box = box #  TODO: Implemented in box should return an array of 3 dimensions
        self._coordinates = coordinates

    @property  # getter
    def coordinates(self):
        return self._particles

    @coordinates.setter #setter
    def num_particles(self):
        self._num_particles = len(self._coordinates)

    @property  # getter
    def coordinates(self):
        return self._particles

    @abstractmethod
    def coordinates(self):
        pass
'''
    @abstractmethod
    def generate_initial_state():
        pass'''

class RandomParticles(ParticlesGenerator):
    ''' This class returns a Particles Object with random positions given number of particles an box dimensions. 
    '''
    def __init__(self, num_particles, box):
        ''' RandomGenerator Class Constructor.

        Parameters:
        -----------
            num_particles (int) : Number of particles of the system to evaluate. 
            box (numpy.array) : An array of shape (3,) with the dimensions of the box (x, y, z). 
        '''
    
    @staticmethod
    def coordinates(self):
        coordinates = (0.5 - np.random.rand(self._num_particles, 3)) * self._box[np.newaxis,:]
        self._coordinates = coordinates


class FileParticles(ParticlesGenerator):
    def __init__(self, file_name):
        self.file_name = file_name
        # self.box = get_box_length() # TODO: Implemented in box

    ''' @Brian
    def get_box_length(self):
        length = self.coordinates.max() - self.coordinates.min()
        return length
        '''
    
    def generate_initial_state(self):
        coordinates = np.loadtxt(self.file_name, skiprows=2, usecols=(1, 2, 3))
        return Particles(coordinates)


class Particles:
    def __init__(self, coordinates):
        self._coordinates = coordinates
       # self._num_particles = num_particles

    @property  # getter
    def coordinates(self):
        return self._coordinates

    @coordinates.setter #setter
    def num_particles(self):
        self._num_particles = len(self._coordinates)

    @property  # getter
    def num_particles(self):
        return self._num_particles

    ''' We really need the user pass coordinates directly to th object?
    @coordinates.setter  # setter
    def coordinates(self, new_coordinates):
        self._coordinates = new_coordinates
        self.bonds = self.build_bond_list()
        '''

'''
    #@property
    def coordinates(self, **kwargs):
        if self.generate_method == 'random':
            return RandomGenerator(self, **kwargs).generate_initial_state()
        else:
            # rise error if file doesn't exist
            return FileGenerator(self, **kwargs).generate_initial_state()
'''

num = 10
box = np.array([10, 10, 10])
part =  RandomGenerator(num, box)
print( part.coordinates )

#print( FileGenerator('./sample_config1.xyz'))
