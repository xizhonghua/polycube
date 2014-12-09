file(GLOB header *.h)
file(GLOB source *.cpp)


link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
$ENV{HOME}/usr/knitro-8.1.1-z/lib
)

message($ENV{HOME}/usr/knitro-8.1.1-z/lib)
#find_package(Boost 1.46 COMPONENTS regex)

add_library(solver SHARED ${header} ${source})
target_link_libraries(solver
${Boost_LIBRARIES}
jtf-optimizer
jtf-function
zjucad-optimizer
ipopt
#blas

#$ENV{HOME}/usr/knitro-8.1.1-z/lib/libknitro.so
#coinmumps
#coinasl
coinhsl
coinmetis
#lapack
#/home/tfjiang/usr/lib/libnlopt.so
#dl
)


