# Range Mutual Information

This package defines functions that compute the amount of information that is expected to be gained about an incomplete occupancy map from a potential range measurement. Range measurement sources could include lidar, sonar, depth cameras, etc. Below is an example output produced by this code.

![mi_surface](https://live.staticflickr.com/65535/49493728457_f33ba30d11_o_d.png)

The package can be built standalone or as a ROS library.

## Building without ROS

Clone the library and build with cmake:

    git clone https://github.com/sportdeath/range_mi.git
    cd range_mi
    mkdir build
    cd build
    cmake ..
    make

## Building with ROS

Clone the library into your catkin workspace:

    cd ~/catkin_ws/src
    git clone https://github.com/sportdeath/range_mi.git
    cd ~/catkin_ws
    catkin_make
