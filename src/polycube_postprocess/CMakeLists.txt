aux_source_directory(. source)
file(GLOB head *.h)

LINK_DIRECTORIES(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_executable(polycube_postprocess ${source} ${head})

target_link_libraries(polycube_postprocess
    foldfree_parameterization
    hj-math
    zjucad-ptree
    zjucad-optimizer
    tetmesh
    jtf-mesh
    jtf-util
    mesh_func
    ${Boost_LIBRARIES}
  )
