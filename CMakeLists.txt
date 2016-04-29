cmake_minimum_required(VERSION 3.4)

enable_language(Fortran)

project(treegrid)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_Fortran_FLAGS "-O -traceback -mcmodel=large -CB")

include_directories("/home/502/aph502/code/fortran/cmake_test/src")
link_directories("/home/502/aph502/code/fortran/cmake_test/src")

# Create an executable file 
add_executable("treegrid" "treegrid.f90" regrid_functions.f90)

target_link_libraries(treegrid modules netcdff)

# Create an executable file 
add_executable("partial_cell_check" "partial_cell_check.f90")

target_link_libraries(partial_cell_check modules netcdff)

# Create an executable file 
add_executable("interpfield" "interpfield.f90" "spline.f90")

target_link_libraries(interpfield modules netcdff)

