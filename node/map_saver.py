#!/usr/bin/env python2

import matplotlib.pyplot as plt
import numpy as np

import rospy
from nav_msgs.msg import OccupancyGrid

class MapSaver:

    def __init__(self):
        topic = "/p_not_measured_map"

        # Subscribe to vectorized maps
        self.sub = rospy.Subscriber(
                topic,
                OccupancyGrid,
                self.map_callback,
                queue_size=1)

    def map_callback(self, map_):
        # Plot the map!
        Z = np.array(map_.data).reshape(map_.info.height, map_.info.width)
        Z = 255 - Z
        plt.imshow(Z, origin='lower', cmap='gray')

        # Get rid of axes
        plt.box(False)
        plt.yticks([])
        plt.xticks([])
        plt.tick_params(direction='in')

        # Save it!
        plt.savefig("occupancy_grid.pdf", bbox_inches='tight', pad_inches=-0.03, transparent=False)

if __name__ == "__main__":
    rospy.init_node("map_saver")
    ms = MapSaver()
    rospy.spin()
