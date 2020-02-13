# Range Mutual Information

This package defines functions that compute the amount of information that is expected to be gained about an incomplete occupancy map from a potential range measurement. Range measurement sources could include lidar, sonar, depth cameras, etc. Below is an example output produced by this code.

![mi_surface](https://live.staticflickr.com/65535/49493728457_f33ba30d11_o_d.png)

The package can be built standalone or as a ROS library.

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

The provided tests compare the accuracy of the provided algorithms to numerical integration as well as providing a timing benchmark. Run them as follows:

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

template <unsigned int dimension, bool lower_bound=true>
void line(
    const unsigned int * const line,
    const double * const vacancy,
    const double * const p_not_measured,
    const double * const width,
    unsigned int num_cells,
    double dtheta,
    double * const output);
