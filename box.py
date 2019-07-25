import numpy as np

class Box:
    def __init__(self, box_dim):
        self.box_dim = box_dim
        
    # wrapping equation
    def wrap(self, box_dim, coordinates):
        """Wraps the coordinates within the box dimensions

        Parameters
        ----------
        box_dim : np.array
            Dimensions of the box in reduced units.
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
        
    def minimum_image_distance(self,coord_i, coord_j, box_dim):
        """Calculate the minimum distance between two atoms.
        
        Parameters
        ----------
        coord_i, coord_j : np.array
            Arrays of the atomic coordinates.
        box_dim : np.array
            Dimensions of the box in reduced units.
        
        Returns
        -------
        coord_ij2 : float
            A scalar product of the positions for two atoms.
        """
        coord_ij = coord_i - coord_j
        coord_ij = coord_ij - box_dim * np.round(coord_ij / box_dim)
        coord_ij2 = np.dot(coord_ij, coord_ij)
        return coord_ij2
