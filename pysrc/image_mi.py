#!/usr/bin/env python3

"""
Convert a folder full of occupancy map
images to MI surface images.
"""

import sys, os
from PIL import Image
import numpy as np
from grid_mi import grid_mi
from matplotlib import cm

# Fetch arguments
if len(sys.argv) != 5:
    print("Usage: ./maps_to_mi num_beams resolution input_dir output_dir")
    sys.exit()

num_beams = int(sys.argv[1])
resolution = float(sys.argv[2])
input_dir = sys.argv[3]
output_dir = sys.argv[4]

# Iterate over files
i = 0
for f in os.listdir(input_dir):
    # Open the array and convert to numpy
    im = Image.open(os.path.join(input_dir, f))
    m = np.array(im)

    # Only take the first 2 dimensions
    while len(m.shape) > 2:
        m = np.take(m, -1, axis=-1)

    # Convert to double
    m = m.astype(np.double)/(2**8 - 1.)

    # Account for resolution
    m = np.power(m, resolution)

    # Compute the MI
    mi = grid_mi(m, num_beams)

    # Normalize
    mi = mi/np.max(mi)

    # Turn back into an image
    mi_im = Image.fromarray((cm.inferno(mi) * (2**8 - 1)).astype(np.uint8))
    mi_im.save(os.path.join(output_dir, f))

    print(i)
    i += 1
