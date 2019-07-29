import numpy as np


class Box:
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
            coordinates -= self.box_dims * \
                np.round(coordinates / self.box_dims[np.newaxis, :])
        return coordinates

    def minimum_image_distance(self, coord_i, coordinates):
        """Calculate the minimum distance between two atoms.

        Parameters
        ----------
        coord_i : np.array
            xyz coordinate of the i-th particle.

        coordinates : np.array
            Array of the atomic xyz coordinate for all particles.

        Returns
        -------
        coord_ij2 : np.array
            Array of the distances between each i-th particle and remaining
            particles
        """
        coord_ij = coord_i[np.newaxis, :] - coordinates[
            coordinates != coord_i].reshape((-1, 3))
        coord_ij = coord_ij - \
            self.box_dims[np.newaxis, :] * \
            np.round(coord_ij / self.box_dims[np.newaxis, :])
        coord_ij2 = np.sum(np.square(coord_ij), axis=1)
        return coord_ij2
