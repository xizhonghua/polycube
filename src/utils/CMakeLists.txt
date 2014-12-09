aux_source_directory(. source)

if(WIN32)
  file(GLOB header_files . "*.h")
  add_executable(hex_utils ${source} ${header_files})
else(WIN32)
  file(GLOB header_files . "*.h")
#include_directories(/opt/localsolver_4_5/include)
link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
#/opt/localsolver_4_5/bin
)

  add_executable(hex_utils ${source})
endif(WIN32)

find_package(Boost 1.46 COMPONENTS filesystem system)
set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

target_link_libraries(hex_utils
  common
  hj_ANN_c
  tetmesh
  hexmesh
  quadmesh
  fast_marching
  hj-math
  zjucad-optimizer
  zjucad-linear-solver
  lapack
  zjucad-ptree
  hex_param
  tet_mesh_sxx
  jtf-mesh
  jtf-util
  ${Boost_LIBRARIES}
  hex_frame_opt
  rot_cubic_f_SH
  #localsolver
  )
