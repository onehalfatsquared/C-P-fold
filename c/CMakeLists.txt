cmake_minimum_required(VERSION 2.8)
project(folding)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp ")

add_subdirectory(/home/anthony/libraries/dlib dlib_build)

add_subdirectory(nauty)
include_directories(nauty)

add_subdirectory(support)
include_directories(support)

add_subdirectory(physics)
include_directories(physics)

add_subdirectory(visual)
include_directories(visual)

add_subdirectory(tpt)
include_directories(tpt)

add_subdirectory(design)
include_directories(design)

add_subdirectory(non_eq_protocol)
include_directories(non_eq_protocol)

add_subdirectory(hydro)
include_directories(hydro)

add_subdirectory(genetic)
include_directories(genetic)

add_subdirectory(lattice)
include_directories(lattice)

add_subdirectory(exe)

file(COPY input DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/exe)
file(COPY dot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/exe)
