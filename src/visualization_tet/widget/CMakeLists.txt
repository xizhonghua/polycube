cmake_minimum_required(VERSION 2.8)
INCLUDE_DIRECTORIES(/home/liquan/usr/include)

#SET(DIFF_SOURCES *.cpp)
#SET(DIFF_HEADER  *.h)
FILE(GLOB HEADFILE *.h)
FILE(GLOB SOURCEFILE *.cpp)
ADD_LIBRARY(lq_widget STATIC ${HEADFILE} ${SOURCEFILE})
TARGET_LINK_LIBRARIES(lq_widget
	vtkWidgets
	vtkHybrid
)
