cmake_minimum_required(VERSION 3.1)

project(range_entropy)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")


######################################
# Compile the library >
######################################

# Add includes
include_directories(include)
file(GLOB SRC_FILES src/*.cpp)
add_library(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
set(LIBS ${LIBS} ${PROJECT_NAME})

# Install the library to CMAKE_INSTALL_PREFIX
# which defaults to /usr/local
install(TARGETS ${PROJECT_NAME}
        DESTINATION lib)
install(DIRECTORY include/${PROJECT_NAME}
        DESTINATION include)

######################################
# < End compile the library
######################################


######################################
# Add tests and ROS >
######################################

option(BUILD_TESTS "BUILD_TESTS" OFF)
if (BUILD_TESTS)
  include_directories(test)
  add_subdirectory(test)
endif()

if(DEFINED CATKIN_DEVEL_PREFIX)
  # Find build dependencies
  find_package(catkin REQUIRED COMPONENTS
    roscpp
    nav_msgs
    geometry_msgs
    message_generation
    )

  # Add the messages
  file(GLOB MSG_FILES msg/*.msg)
  foreach(_file ${MSG_FILES})
    get_filename_component(_name ${_file} NAME)
    add_message_files(
      FILES
      ${_name}
    )
  endforeach()

  # Generate messages
  generate_messages(
    DEPENDENCIES
    std_msgs
  )

  # Define the package
  catkin_package(
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME}
    CATKIN_DEPENDS message_runtime
  )

  # Add nodes
  include_directories(node)
  add_subdirectory(node)

  # Install the library
  install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )
  install(DIRECTORY include/${PROJECT_NAME}/
    DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  )
endif()

######################################
# < End add tests and ROS
######################################
