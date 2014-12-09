cmake_minimum_required(VERSION 2.8)
INCLUDE_DIRECTORIES(/home/liquan/usr/include)	
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
#INCLUDE_DIRECTORIES(../tetmesh)
FILE(GLOB SimpleViewHeaders *.h)
FILE(GLOB SimpleViewUI  *.ui)
FILE(GLOB SimpleViewSrcs  *.cpp)

QT4_WRAP_UI(UISrcs ${SimpleViewUI})
QT4_WRAP_CPP(MOCSrcs  ${SimpleViewHeaders})

ADD_LIBRARY(main_win ${SimpleViewSrcs} ${UISrcs} ${MOCSrcs})
TARGET_LINK_LIBRARIES(main_win
	${VTK_LIBRARIES} 
	${QT_LIBRARIES}
	side_bar
	call_back
	lq_widget
	io
	vtk_operate
	config
	QtCore
	QVTK
	tetmesh
	common
	lapack)
