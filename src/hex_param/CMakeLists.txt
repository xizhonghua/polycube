file(GLOB headers *.h)
aux_source_directory(. source)

add_library(hex_param STATIC ${source} ${headers})

set(Boost_USE_STATIC_LIBS        OFF)  # I've already tried ON
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(Boost 1.46 COMPONENTS filesystem system)
message("boost_libraries ${Boost_LIBRARIES}")
set(CMAKE_CXX_FLAGS "-std=c++0x  -fPIC ${CMAKE_CXX_FLAGS}")
target_link_libraries(hex_param
  common
  numeric
  hj-function
  jtf-util
  equation_graph
  ${Boost_LIBRAIES}
  )
