#!/usr/bin/env python3

from PIL import Image
import numpy as np

occupied_threshold = 0.5
shrink_factor = 3
image_file = "building_31.png"

m = np.array(Image.open(image_file))

new_height = int(m.shape[0]/shrink_factor)
new_width  = int(m.shape[1]/shrink_factor)

m_shrunk = np.ones((new_height, new_width))

for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        vacancy = m[i,j]/255.
        if vacancy < occupied_threshold:
            m_shrunk[int(i/shrink_factor),int(j/shrink_factor)] = 0

m_shrunk = Image.fromarray(np.uint8(m_shrunk*255))
m_shrunk.save("shrunk_" + image_file)
