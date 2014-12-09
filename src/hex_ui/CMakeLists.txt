aux_source_directory(. source)
file(GLOB headers *.h)

link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_executable(hex_ui ${source} ${headers})

find_package(Boost 1.46 COMPONENTS filesystem system)
#message("${Boost_LIBRARIES}")
#add_definitions(-fopenmp -fpermissive -fPIC)
#link_directories(/usr/lib)


set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
target_link_libraries(hex_ui
  hex_frame_opt
  zjucad-ptree
  common
  tetmesh
  tet_mesh_sxx
  hexmesh
  trimesh
  zjucad-optimizer
  minpack
  boost_thread
  rot_cubic_f_SH
  hex_param
  hj_ANN_c
  hex_process
  tetmesh_refine
  hex_process
  jtf-mesh
  jtf-util
  ${Boost_LIBRARIES}
  ipopt
  coinhsl
  coinmetis
  jtf-function
  fairing
  mesh_func
  )
