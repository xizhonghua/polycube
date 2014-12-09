aux_source_directory(. source)
file(GLOB header *.h)

add_library(numeric SHARED ${source} ${header})

target_link_libraries(numeric
zjucad-optimizer
lapack
${Boost_LIBRARIES})
