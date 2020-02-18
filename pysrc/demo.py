#!/usr/bin/env python3

import numpy as np
from grid_mi import grid_mi

# Create an occupancy map
map_ = np.array(
       [[0.8,  0.8,  0.8,  0.8,  0.8,  0.8, 0.8],
        [0.8,  0.8,  0.8,  0.8,  0.8,  0.8, 0.8],
        [0.8,  0.8,  0.8,  0.8,  0.8,  0.8, 0.8],
        [0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8],
        [0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8],
        [0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8],
        [0.8, 0.01, 0.99, 0.99, 0.99, 0.01, 0.8]])

num_beams = 200

# Compute!
print(grid_mi(map_, num_beams))
