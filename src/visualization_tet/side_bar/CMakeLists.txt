cmake_minimum_required(VERSION 2.8)
FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )
INCLUDE_DIRECTORIES(/home/liquan/usr/include)

FILE(GLOB ViewSrcs *.cpp)
FILE(GLOB ViewUI *.ui)
FILE(GLOB ViewHeaders *.h)

QT4_WRAP_UI(UISrcs ${ViewUI})
QT4_WRAP_CPP(MOCSrcs ${ViewHeaders})

ADD_LIBRARY(side_bar  STATIC ${ViewSrcs} ${UISrcs} ${MOCSrcs})

target_link_libraries(side_bar
  QVTK
  vtkWidgets 
  vtkHybrid
)
