file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(spinxform ${source} ${header})
