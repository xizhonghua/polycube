file(GLOB head *.h)
file(GLOB src *.cpp)
add_library(tetmesh_refine SHARED ${head} ${src})

target_link_libraries(tetmesh_refine
  hj-math
  zjucad-optimizer)
