aux_source_directory(. source)
file(GLOB header *.h)
add_library(hex_process SHARED ${source} ${header})
target_link_libraries(hex_process
    hexmesh
    quadmesh
)
