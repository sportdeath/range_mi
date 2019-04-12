cmake_minimum_required(VERSION 3.1)

project(wandering_robot)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

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
# Compile the tests >
######################################


option(BUILD_TESTS "BUILD_TESTS" OFF)
if (BUILD_TESTS)
  add_subdirectory(test)
endif()

######################################
# < End compile the tests
######################################
