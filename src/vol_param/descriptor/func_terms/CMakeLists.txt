file(GLOB  header_files *.h)
file(GLOB source *.cpp)

add_library(func_terms SHARED  ${header_files} ${source})
target_link_libraries(func_terms
mesh_func)
