#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np

import rospy
from range_mi.msg import MIGrid

class MIGridHeatmap:

    def __init__(self):
        topic = "/mi"

        # Subscribe to vectorized maps
        self.sub = rospy.Subscriber(
                topic,
                MIGrid,
                self.map_callback,
                queue_size=1)

    def map_callback(self, map_):
        # Plot the map!
        Z = np.array(map_.data).reshape(map_.height, map_.width)
        plt.imshow(Z, origin='lower', cmap='inferno')

        # Get rid of axes
        plt.box(False)
        plt.yticks([])
        plt.xticks([])
        plt.tick_params(direction='in')

        # Save it!
        plt.savefig("mi_grid.pdf", bbox_inches='tight', pad_inches=-0.03, transparent=False)
        plt.colorbar()
        plt.savefig("mi_grid_colorbar.pdf", bbox_inches='tight', pad_inches=0, transparent=False)

if __name__ == "__main__":
    rospy.init_node("mi_grid_heatmap")
    hm = MIGridHeatmap()
    rospy.spin()
