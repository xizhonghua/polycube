if(COMMAND cmake_minimum_required)
   cmake_minimum_required(VERSION 2.8)
endif(COMMAND cmake_minimum_required)

PROJECT(foldfree_parameterization)

file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.h")

INCLUDE_DIRECTORIES(
$ENV{HOME}/usr/include
)

LINK_DIRECTORIES(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

#ADD_EXECUTABLE(foldfree_parameterization ${SOURCES} ${HEADERS})
add_library(foldfree_parameterization ${SOURCES} ${HEADERS})

SET(CMAKE_CXX_FLAGS -fopenmp)
SET(CMAKE_CXX_LINKFLAGS -Wl,--rpath -Wl)
set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
TARGET_LINK_LIBRARIES(foldfree_parameterization
zjucad-ptree
ipopt
coinhsl
coinmetis
dl
blas
lapack
${Boost_LIBRARIES})



