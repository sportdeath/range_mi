# Range Mutual Information

This package defines functions that compute the amount of information that is expected to be gained about an incomplete occupancy map from a potential range measurement. Range measurement sources could include lidar, sonar, depth cameras, etc. Below is an example output produced by this code.
The package can be built standalone or as a ROS library.

![mi_surface](https://live.staticflickr.com/65535/49493728457_f33ba30d11_o_d.png)


## Usage

### Without ROS

Clone the library and build with ```cmake```:

    git clone https://github.com/sportdeath/range_mi.git
    cd range_mi
    mkdir build
    cd build
    cmake ..
    make

To build the tests you can add a flag to ```cmake``` before making:

    cmake -DBUILD_TESTS=ON ..
    make

The tests can then be run as follows:

    ./test/api_demo
    ./test/profile_mi
    ./test/numerical_mi
    ./test/numerical_mi_distorted

### With ROS

Clone the library into your catkin workspace:

    cd ~/catkin_ws/src
    git clone https://github.com/sportdeath/range_mi.git
    cd ~/catkin_ws
    catkin_make

## API

A basic demonstration of how this code can be used for computing 2D mutual information surfaces is found in ```test/api_demi.cpp```. Run it as described above in the "Without ROS" section.
