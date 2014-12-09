aux_source_directory(. source)
#INCLUDE_DIRECTORIES(/home/liquan/usr/include)

LINK_DIRECTORIES(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

# to include mesquite
#include_directories(/usr/local/include)
add_library(tetmesh SHARED ${source})
target_link_libraries(
  tetmesh
  jtf-mesh
  hj-math
  jtf-function
  jtf-optimizer
 # mesquite
  common
ipopt
coinhsl
coinmetis
)

