import numpy as np

class Box:
    """Holds all the information for the Box.

    Parameters
    ----------
    box_dim : np.array
        The dimensional lengths of the box, should be a numpy array ([x, y, z]).
        With shape (1, 3).

    Returns
    -------
    self : Box
        Returns an instance of itself.
    
    Attributes
    ----------
    box_dim : np.array
        The dimensional lengths of the box, should be a numpy array ([x, y, z]).
    """
    def __init__(self, box_dim):
        self.box_dim = box_dim

    @property
    def volume(self):
        """Calculate the box volume

        Returns
        -------
        box volume : float
            Computed box volume
        """
        return np.prod(self.box_dim)
        
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
        coord_ij2 = np.sum(np.square(coord_ij), axis = 1)
        return coord_ij2
