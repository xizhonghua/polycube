cmake_minimum_required(VERSION 2.8)
INCLUDE_DIRECTORIES(/home/liquan/usr/include)	
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
SET(main_srcs main.cpp)
set(CMAKE_DEBUG_POSTFIX "d")
link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

ADD_EXECUTABLE(visual_tet ${main_srcs})

TARGET_LINK_LIBRARIES(visual_tet 
	side_bar
	call_back
	lq_widget
	io
	vtk_operate
	config
	main_win
	QtCore
)
