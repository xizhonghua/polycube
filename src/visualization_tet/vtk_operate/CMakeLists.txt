cmake_minimum_required(VERSION 2.8)
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
INCLUDE_DIRECTORIES(/home/liquan/usr/include)
#SET(DIFF_SOURCES *.cpp)
#SET(DIFF_HEADER  *.h)
FILE(GLOB HEADFILE *.h)
FILE(GLOB SOURCEFILE *.cpp)

ADD_LIBRARY(vtk_operate STATIC ${HEADFILE} ${SOURCEFILE})

TARGET_LINK_LIBRARIES(vtk_operate
	call_back
	io
	lq_widget	
	QVTK
	vtkWidgets
    vtkHybrid
	tetmesh
	common
	lapack
)
