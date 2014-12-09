aux_source_directory(. source)
file(GLOB header *.h)
add_library(quadmesh SHARED ${source} ${header})
