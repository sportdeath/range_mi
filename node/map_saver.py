#!/usr/bin/env python2

import numpy as np
from PIL import Image

import rospy
from nav_msgs.msg import OccupancyGrid

class MapSaver:

    def __init__(self):
        topic = "/mapped"

        # Subscribe to vectorized maps
        self.sub = rospy.Subscriber(
                topic,
                OccupancyGrid,
                self.map_callback,
                queue_size=1)

    def map_callback(self, map_):
        # Plot the map!
        Z = np.array(map_.data).reshape(map_.info.height, map_.info.width)
        Z = (1 - Z/100.)
        #Z[np.logical_and(Z < 1, Z > 0)] = 0.5
        Z = np.flip(Z,0)
        Z *= 255
        Z = Z.astype(np.uint8)

        img = Image.fromarray(Z)
        img.save('occupancy_grid.png')
        print("saved!")

if __name__ == "__main__":
    rospy.init_node("map_saver")
    ms = MapSaver()
    rospy.spin()
