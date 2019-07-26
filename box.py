import numpy as np

class Box:
    def __init__(self, box_dim):
        self.box_dim = box_dim
        self.box_vol = box_vol

    def volume(self):
        """Calculate the box volume

        Returns
        -------
        box_vol : float
            Computed box volume
        """
        self.box_vol = np.dot(self.box_dim,self.box_dim)
        return self.box_vol
        
    def wrap(self, coordinates):
        """Wraps the coordinates within the box dimensions

        Parameters
        ----------
        coordinates : np.array
            Array of the atomic coordinates.

        Returns
        -------
        coordinates : np.array
            Arrays of the wrapped atomic coordinates.
        """
        if len(coordinates.shape) == 1:
            coordinates -= self.box_dim * np.round(coordinates / self.box_dim)
        else:
            coordinates -= self.box_dim * np.round(coordinates / self.box_dim[np.newaxis, :])
        return coordinates
        
    def minimum_image_distance(self,coord_i, coordinates):
        """Calculate the minimum distance between two atoms.
        
        Parameters
        ----------
        coord_i : np.array
            xyz coordinate of the i-th particle.

        coord_j : np.array
            Array of the atomic xyz coordinate for all particles.
        
        Returns
        -------
        coord_ij2 : np.array
            Array of the distances between each i-th particle and remaining particles
        """
        coord_ij = coord_i[np.newaxis, :] - coord_j[coord_j != coord_i]
        coord_ij = coord_ij - self.box_dim[np.newaxis, :] * np.round(coord_ij / self.box_dim[np.newaxis, :])
        coord_ij2 = np.square(coord_ij)
        return coord_ij2
