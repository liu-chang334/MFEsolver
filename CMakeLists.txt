# Minimum required CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and settings
project(FEsolver LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Check and set VTK_DIR using VTK_ROOT_FOLDER
if(DEFINED ENV{VTK_ROOT_FOLDER})
    set(VTK_DIR $ENV{VTK_ROOT_FOLDER}/lib/cmake/vtk-9.4)
    message(STATUS "VTK root directory: $ENV{VTK_ROOT}")
    message(STATUS "VTK config directory: ${VTK_DIR}")
else()
    message(FATAL_ERROR "VTK_ROOT environment variable is not defined.")
endif()

# Include Eigen library using an environment variable
if(DEFINED ENV{EIGEN3_ROOT})
    include_directories($ENV{EIGEN3_ROOT})
    message(STATUS "Eigen3 root directory: $ENV{EIGEN3_ROOT}")
else()
    message(FATAL_ERROR "EIGEN3_ROOT environment variable is not defined.")
endif()

# Find the VTK library
find_package(VTK REQUIRED)

if(VTK_FOUND)
    message(STATUS "VTK found: ${VTK_DIR}")
    include(${VTK_USE_FILE})
else()
    message(FATAL_ERROR "VTK not found. Check VTK_DIR or installation.")
endif()

# Add source files
file(GLOB SRC_FILES src/*.cpp)
# file(GLOB SRC_FILES vtksample/vtkPoints.cpp)
# file(GLOB SRC_FILES vtksample/vtkUnstructuredGrid.cpp)
# file(GLOB SRC_FILES vtksample/plotcase3.cpp)

# Add the executable
add_executable(FEsolver main.cpp ${SRC_FILES})
# add_executable(FEsolver ${SRC_FILES})

# Link VTK libraries
target_link_libraries(FEsolver ${VTK_LIBRARIES})

# Enable Unity Build to speed up compilation
set(CMAKE_UNITY_BUILD ON)

# Print compiler information
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
