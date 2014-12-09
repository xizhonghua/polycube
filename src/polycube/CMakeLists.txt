aux_source_directory(. source)
file(GLOB head *.h)

LINK_DIRECTORIES(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_executable(polycube ${source} ${head})

target_link_libraries(polycube
    hj-math
    zjucad-ptree
    zjucad-optimizer
    tetmesh
    jtf-mesh
    tet_mesh_sxx
    jtf-util
    jtf-optimizer
    jtf-function
    mesh_func
    ${Boost_LIBRARIES}
  )
