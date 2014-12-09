cmake_minimum_required(VERSION 2.8)

project(hexgen)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{HOME}/usr/share/cmake/Modules/")
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)

set(CMAKE_CXX_FLAGS "-std=c++0x -fopenmp ${CMAKE_CXX_FLAGS}")
add_subdirectory(src)
