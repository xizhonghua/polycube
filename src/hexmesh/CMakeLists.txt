aux_source_directory(. source)
file(GLOB header *.h)
add_library(hexmesh SHARED ${source} ${header})
target_link_libraries(
  hexmesh
  mesquite
)
