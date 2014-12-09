file(GLOB header *.h)
file(GLOB source *.cpp)

LINK_DIRECTORIES(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_executable(param ${source} ${header})
target_link_libraries(param
solver
descriptor
func_terms
zjucad-ptree
jtf-mesh
tet_mesh_sxx
common
equation_graph
hex_param
param_common
jtf-function
mesh_func
ipopt
coinhsl
coinmetis
)
