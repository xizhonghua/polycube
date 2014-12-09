file(GLOB  header_files *.h)
file(GLOB  source_files *.cpp)

link_directories(/home/tfjiang/usr/Linux/64/gcc4/lib)
add_library(param_common SHARED ${header_files} ${source_files})
target_link_libraries(param_common jtf-mesh)

#add_custom_target(param_common SOURCES ${header_files})
