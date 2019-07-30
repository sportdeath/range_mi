#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np

import rospy
from range_entropy.msg import EntropyMap

class Heatmap:

    def __init__(self):
        topic = "/entropy"

        # Subscribe to vectorized maps
        self.sub = rospy.Subscriber(
                topic,
                EntropyMap,
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
        plt.savefig("entropy_surface.pdf", bbox_inches='tight', pad_inches=-0.03, transparent=False)
        plt.colorbar()
        plt.savefig("entropy_surface_colorbar.pdf", bbox_inches='tight', pad_inches=0, transparent=False)

if __name__ == "__main__":
    rospy.init_node("heatmap")
    hm = Heatmap()
    rospy.spin()
