import numpy as np

class RandomParticles:
    ''' This class returns a Particles Object with random positions given number of particles an box dimensions. 
    '''
    def __init__(self, num_particles, box):
        #super().__init__()
        
        self._num_particles = num_particles
        self._box = box
        self._coordinates = self.get_coordinates()
        
        ''' RandomGenerator Class Constructor.

        Parameters:
        -----------
            num_particles (int) : Number of particles of the system to evaluate. 
            box (numpy.array) : An array of shape (3,) with the dimensions of the box (x, y, z). 
        '''

    @property
    def coordinates(self):
        return self._coordinates
    
    # 
    def get_coordinates(self):
        a = (0.5 - np.random.rand(self._num_particles, 3)) * self._box[np.newaxis,:]
        return a

    def num_particles(self):
        return self._num_particles

    def get_num_particles(self):
        self._num_particles = len(self._coordinates)

num = 10
box = np.array([10, 10, 10])
part =  RandomParticles(num, box)
print( part.coordinates )