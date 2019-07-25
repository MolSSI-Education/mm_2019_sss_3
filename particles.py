import numpy as np
from abc import ABC, abstractmethod

# Interfase
class ParticlesGenerator(ABC):
    @abstractmethod
    def generate_initial_state():
        pass

class RandomGenerator(ParticlesGenerator):
    def __init__(self, Particles, num_particles, box_length):
        # assuming the object Box has an atribute or a
        self.num_particles = num_particles
        self.box_length = box_length

    def generate_initial_state(self, Particles):
        coordinates = (0.5 - np.random.rand(self.num_particles)) * self.box_length
        return coordinates

class FileGenerator(ParticlesGenerator):
    def __init__(self, Particles, file_name):
        self.file_name = file_name
        # self.box_length # Implemented in box

    @property
    def box_length(self):
        length = self.coordinates.max() - self.coordinates.min()
        return length
    
    def generate_initial_state(self, Particles):
        Particles.coordinates = np.loadtxt(self.file_name, skiprows=2, usecols=(1, 2, 3))
''' 
Brian
    def get_box_length(self):
        box_length = 
'''

class Particles:
    def __init__(self, generate_method):
        self.generate_method = generate_method
        self.coordinates

    @property
    def coordinates(self, **kwargs):
        if self.generate_method == 'random':
            return RandomGenerator(self, **kwargs).generate_initial_state()
        else:
            # rise error if file doesn't exist
            return FileGenerator(self, **kwargs).generate_initial_state()


