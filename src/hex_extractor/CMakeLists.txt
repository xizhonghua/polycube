file(GLOB header *.h)
file(GLOB source *.cpp)

find_package(Boost 1.46 COMPONENTS filesystem system)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
else()
    message("can not find boost")
endif()


add_executable(test_extract_integer_points ${header} ${source})
target_link_libraries(test_extract_integer_points
lapack
jtf-mesh
jtf-util
${Boost_LIBRARIES}
)

