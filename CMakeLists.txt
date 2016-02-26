cmake_minimum_required(VERSION 3.4)

enable_language(Fortran)

project(treegrid)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories("/home/502/aph502/code/fortran/cmake_test/src")
link_directories("/home/502/aph502/code/fortran/cmake_test/src")

# Create an executable file called helloworld from sources:
add_executable("treegrid" "treegrid.f90")

target_link_libraries(treegrid modules netcdff)
