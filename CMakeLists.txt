# Minimum CMake version & Project
cmake_minimum_required(VERSION 3.16)
project(MFEsolver LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Generator platform
if(MSVC AND CMAKE_GENERATOR_PLATFORM)
  message(STATUS "Generator platform: ${CMAKE_GENERATOR_PLATFORM}")
endif()

# src
file(GLOB_RECURSE ELEM_SRC    "elem/*.cpp")
file(GLOB_RECURSE SOLVER_SRC  "solver/*.cpp")
file(GLOB_RECURSE UTILS_SRC   "utils/*.cpp")
file(GLOB_RECURSE CONS_SRC     "cons/*.cpp")

set(SRC_FILES
    main.cpp
    ${ELEM_SRC}
    ${SOLVER_SRC}
    ${UTILS_SRC}
    ${CONS_SRC}
)

# executable
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}       # root dir
    ${CMAKE_SOURCE_DIR}/elem  # elem 
    ${CMAKE_SOURCE_DIR}/cons  # cons

)

# Eigen3
target_include_directories(${PROJECT_NAME} PRIVATE $ENV{EIGEN3_ROOT})