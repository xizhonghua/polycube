project(fairing)
file(GLOB header *.h)
file(GLOB source *.cpp)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Werror -ansi -pedantic")
#add_executable(fairing ${header} ${source})
add_library(fairing STATIC ${source} ${header})

target_link_libraries(fairing 
spqr
lapack
/usr/lib/libblas.so.3
cholmod 
umfpack
colamd
ccolamd 
camd 
amd 
m
)
