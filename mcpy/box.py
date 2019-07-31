import numpy as np


class Box:
    """Holds all the information for the Box.

    Parameters
    ----------
    box_dims : np.array
        The dimensional lengths of the box, should be a numpy array ([x, y, z]).
        With shape (1, 3).

    Returns
    -------
    self : Box
        Returns an instance of itself.
    
    Attributes
    ----------
    box_dims : np.array
        The dimensional lengths of the box, should be a numpy array ([x, y, z]).
    """
    def __init__(self, box_dims):
        self.box_dims = box_dims

    @property
    def volume(self):
        """Calculate the box volume

        Returns
        -------
        box volume : float
            Computed box volume
        """
        return np.prod(self.box_dims)

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
            coordinates -= self.box_dims * \
                np.round(coordinates / self.box_dims)
        else:
            coordinates -= self.box_dims[np.newaxis, :] * \
                np.round(coordinates / self.box_dims[np.newaxis, :])
        return coordinates

    def minimum_image_distance(self, index, coordinates):
        """Calculate the minimum distance between two atoms.

        Parameters
        ----------
        index : int
            index of the particle to take the minimum images for

        coordinates : np.array
            Array of the atomic xyz coordinate for all particles.

        Returns
        -------
        coord_ij2 : np.array
            Array of the distances between each i-th particle and remaining
            particles
        """
        if index != 0 and index != len(coordinates):
            coord_ij = coordinates[index, :] - coordinates[:index, :]
            temp = coordinates[index, :] - coordinates[index + 1:, :]
            coord_ij = np.concatenate((coord_ij, temp))
        elif index == len(coordinates):
            coord_ij = coordinates[index, :] - coordinates[:index, :]
        elif index == 0:
            coord_ij = coordinates[index, :] - coordinates[index + 1:, :]
        
        coord_ij = coord_ij - \
            self.box_dims[np.newaxis, :] * \
            np.round(coord_ij / self.box_dims[np.newaxis, :])
        coord_ij2 = np.sum(np.square(coord_ij), axis=1)
        return coord_ij2
