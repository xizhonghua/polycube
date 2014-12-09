aux_source_directory(. source)

link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_executable(deformer ${source})

target_link_libraries(deformer
	zjucad-ptree
  zjucad-optimizer
  zjucad-linear-solver
        #tet-vol
	tetmesh
	trimesh
	common
  lapack
  ${Boost_LIBRARIES}
  )
