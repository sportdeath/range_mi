# Include catkin_libraries
include_directories(${catkin_INCLUDE_DIRS})
set(LIBS ${LIBS} ${catkin_LIBRARIES})

# Add the nodes
file(GLOB SRC_FILES *.cpp)
foreach(_file ${SRC_FILES})
  get_filename_component(_name ${_file} NAME_WE)
  add_executable(${_name} ${_file})
  add_dependencies(${_name} ${${PROJECT_NAME}_EXPORTED_TARGETS})
  target_link_libraries(${_name} ${LIBS})
endforeach()
