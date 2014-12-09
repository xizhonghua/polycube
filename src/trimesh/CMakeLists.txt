aux_source_directory(. source)
file(GLOB header *.h)
add_library(trimesh SHARED ${source} ${header})
