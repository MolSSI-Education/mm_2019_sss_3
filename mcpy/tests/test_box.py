from mcpy import Box
import pytest
import sys
import numpy as np


@pytest.mark.parametrize("box_dims, volume", [([2, 2, 2], 8),
                                              ([1, 2, 3], 6),
                                              ([1, 10, 0], 0)])
def test_volume_calculation(box_dims, volume):
    calculated_volume = Box.volume(box_dims)
    assert np.isclose(calculated_volume, volume), \
        "The calculated value does not match the box volume"


def






