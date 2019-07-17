#!/usr/bin/env python2

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
        x = np.arange(map_.width)
        y = np.arange(map_.height)
        X, Y = np.meshgrid(x, y)

        Z = np.array(map_.data).reshape(map_.width, map_.height)

        fig = plt.figure()
        a = fig.add_subplot(111, projection='3d')
        a.plot_surface(X, Y, Z, cmap='inferno')
        plt.show()

if __name__ == "__main__":
    rospy.init_node("heatmap")
    hm = Heatmap()
    rospy.spin()
