"""
particles.py
A python package for the MolSSI Software Summer School.
Contains Particles class
"""

import numpy as np
from box import *

class Particles():
    def __init__(self, file_name = None, num_particles=0, box_dims=None):
        ''' Particles Class Constructor.

        Parameters:
        -----------
            file_name : str
                A string with the path to the file with x, y, z coordinates.
            num_particles : int
                Number of particles of the system to generate ramdomly if file_name is not given.
            box_dims : np.array
                Array of shape (3,) of x, y, z dimensions.
        '''
        if file_name is None:
            if num_particles == 0 or box_dims is None :
                # TODO: Edit
                print('Please provide a file_name with particles coordinates, or\nthe number of particles and box dimensions.')
            else:
            # Generate particles ramdomly
                self._num_particles = num_particles
                self._box = Box(box_dims)
                self._coordinates = (0.5 - np.random.rand(self._num_particles, 3)) * self._box.box_dim[np.newaxis,:]
        else:
            # Generate particles from file
            self._coordinates = np.loadtxt(file_name, skiprows=2, usecols=(1, 2, 3))
            self._num_particles = len(self._coordinates)
            self._box = self._create_box(box_dims)

    def __str__(self):
        return( F'Particles Object: {self._num_particles} particles.' )

    def _create_box(self, box_dims):
        '''Generates the Box object the FetchParticles Object

            If box_dims is an np.array the coordinates are wrapped accordingly to the box_dims.
            Else, if box_dims is None, the box object dimensions are taken from the coordinates given as a file.

            Parameters:
            -----------
            box_dims : np.array, None
                Array of shape (3,) of x, y, z dimensions.
        
            Returns:
            --------
            box : Box class object
                An Box class object with given dimensions.'''
        if box_dims is None:
            box_dims = abs( self.coordinates.max(axis=0) - self.coordinates.min(axis=0) )
            box = Box(box_dims)

        elif isinstance(box_dims, np.ndarray):
            box = Box(box_dims)
            self._coordinates = box.wrap(self._coordinates)

        return box

    @property
    def coordinates(self):
        '''Returns the coordinates of the Particles Object'''
        return self._coordinates

    @coordinates.setter
    def coordinates(self, coordinates):
        '''Sets the coordinates of the Particles Object'''
        self._coordinates = coordinates

    @property
    def num_particles(self):
        '''Returns the coordinates of the Particles Object'''
        return self._num_particles

    @num_particles.setter
    def num_particles(self, num_particles):
        '''Sets the coordinates of the Particles Object'''
        self._num_particles = num_particles

    @property
    def box(self):
        '''Returns the Box Object of the Particles Object'''
        return self._box
 
'''
num_particles = 3
box_dims = np.array([1, 1, 1])

# If we load from file and want to wrap
part1 =  Particles(file_name='sample_config1.xyz', box_dims= box_dims)
# If we load from file with given dimensions
part2 =  Particles(file_name='sample_config1.xyz')
# If we create them randomly
part3 =  Particles( num_particles = 3, box_dims = box_dims)

print( part3.coordinates)
print(F'\nNum particles: {part3.num_particles} \nBox dims: {part3.box.box_dim}')
'''