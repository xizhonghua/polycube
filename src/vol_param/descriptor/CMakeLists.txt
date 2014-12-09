file(GLOB header *.h)
file(GLOB source *.cpp)

add_subdirectory(func_terms)
link_directories($ENV{HOME}/usr/Linux/64/gcc4/lib)
add_library(descriptor SHARED ${header} ${source})
target_link_libraries(descriptor
func_terms
hj-math
#equation_graph
mesh_func
jtf-function
)
