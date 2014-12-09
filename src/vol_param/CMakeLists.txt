set(CMAKE_CXX_FLAGS "-fopenmp -fPIC")

include_geo_sim_sdk()
link_geo_sim_sdk()

set(CMAKE_DEBUG_POSTFIX "d")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.46 COMPONENTS filesystem system serialization date_time)
set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

add_subdirectory(common)
add_subdirectory(descriptor)
add_subdirectory(param)
add_subdirectory(solver)
