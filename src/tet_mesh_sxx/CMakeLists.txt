aux_source_directory(. source)
file(GLOB header *.h)

add_library(tet_mesh_sxx SHARED ${source} ${header})
