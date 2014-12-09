#project(test)
set(MY_MODULE_PATH $ENV{HOME}/usr/share/cmake/Modules)

list(APPEND CMAKE_MODULE_PATH ${MY_MODULE_PATH})

find_package(CppUnit)

if(CPPUNIT_FOUND)
  include_directories(${CPPUNIT_INCLUDE_DIR})
else()
  message("can not find cppunit")
endif()

#set(Boost_USE_STATIC_LIBS        OFF)
#set(Boost_USE_MULTITHREADED      ON)
#set(Boost_USE_STATIC_RUNTIME    OFF)

#find_package( Boost 1.42 COMPONENTS date_time system program_options )

#if(Boost_FOUND)
##  include_directories(${Boost_INCLUDE_DIRS})
#  message(${Boost_LIBRARIES})
#else()
#  message("--[info] can not find Boost")
#endif()


file(GLOB source *.cpp)
file(GLOB header *.h)

add_executable(test ${source} ${header})

target_link_libraries(test
  ${CPPUNIT_LIBRARIES}
  hex_param
  tetmesh
  hexmesh
  common
  trimesh
  numeric
  zjucad-ptree
  ${Boost_LIBRARIES}
)
