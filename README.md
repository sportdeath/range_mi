# Range Mutual Information

This package defines functions that compute the amount of information that is expected to be gained about an incomplete occupancy map from a potential range measurement. Range measurement sources could include lidar, sonar, depth cameras, etc. Below is an example output produced by this code.
The package can be built standalone or as a ROS library.

![mi_surface](https://live.staticflickr.com/65535/49493728457_f33ba30d11_o_d.png)


## Usage

### C++ API

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

To run a test, for example the ```test/api_demo.cpp``` test, execute:

    ./test/api_demo

Taking a look at the ```test/api_demo.cpp``` script is a good place to start to understand the API for computing 2D mutual information surfaces.

### Python API

After building the C++ library as described above also install the library:

    sudo make install

Then compile the python library (requires ```cython```):

    python setup.py build_ext --inplace

This builds a function ```grid_mi``` that simply takes as input a 2D occupancy map and returns a 2D mutual information surface (both numpy arrays).

### ROS Demo Nodes

Clone the library into your catkin workspace:

    cd ~/catkin_ws/src
    git clone https://github.com/sportdeath/range_mi.git
    cd ~/catkin_ws
    catkin_make

Then to compute the mutual information of an occupancy grid, you can run:

    roslaunch range_mi occupancy_grid_mi.launch

In ```rviz``` the occupancy grid map will be displayed on the topic ```/map``` and the mutual information surface will be displayed on the topic ```/mi_map```. The topics names, number of range measurement beams and other parameters are available in the ```launch/params.yaml``` file. By running the ```node/mi_grid_heatmap.py``` script, the mutual information surface can be saved as colored image like the one above.

By using the "Publish Point" tooltip in ```rviz```, you can click on a point in the map to condition the mutual information on a measurement being made at that point. This will likely reduce the information gain within the view of the clicked point.

We have also included another script that computes the mutual information by averaging the information gain over random occupancy assigments. This will take much longer than the included algorithm, but it will eventually converge to the same result. Run it on a map as follows:

    roslaunch range_mi monte_carlo_mi.launch
