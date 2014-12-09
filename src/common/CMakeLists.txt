#aux_source_directory(. source)

INCLUDE_DIRECTORIES(/home/liquan/usr/include)
file(GLOB source *.cpp)
file(GLOB header_files *.h)

link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

if(WIN32)
  add_definitions(-DHEXGEN_COMMON_EXPORT)
  file(GLOB header_files . "*.h")
  add_library(common SHARED ${source} ${header_files})
else(WIN32)
  add_library(common SHARED ${source} ${header_files})
endif(WIN32)

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
#add_subdirectory(hjlib)

target_link_libraries(common
  hj-function
  visit_verdict
  numeric
  ipopt
  coinhsl
  coinmetis
)
