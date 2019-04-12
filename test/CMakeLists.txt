# Enable testing
enable_testing()

# Compile all the test binaries
file(GLOB SRC_FILES *.cpp)
foreach(_file ${SRC_FILES})
    get_filename_component(_name ${_file} NAME_WE)
    add_executable(${_name} ${_file})
    target_link_libraries(${_name} ${LIBS})
    add_test(${_name} ${_name})
endforeach()
