cmake_minimum_required(VERSION 2.8)
INCLUDE_DIRECTORIES(/home/liquan/usr/include)

#SET(DIFF_SOURCES *.cpp)
#SET(DIFF_HEADER  *.h)
FILE(GLOB HEADFILE *.h)
FILE(GLOB SOURCEFILE *.cpp)
ADD_LIBRARY(call_back STATIC ${HEADFILE} ${SOURCEFILE})
TARGET_LINK_LIBRARIES(call_back
	vtkWidgets
	vtkHybrid
)
