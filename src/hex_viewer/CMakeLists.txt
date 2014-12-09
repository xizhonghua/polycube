if(WIN32)
  add_definitions(-D_WCHAR_T_DEFINED)
  set(GLUT_LIBRARIES glut32)
  set(OPENSCENEGRAPH_LIBRARIES osg osgViewer osgGA)
else(WIN32)
  find_package(OpenSceneGraph REQUIRED osg osgViewer osgGA)
  find_package(GLUT)
endif(WIN32)

aux_source_directory(. source)
include_directories(/usr/include/eigen3)
include_directories($ENV{HOME}/usr/include/AntTweakBar)
link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_executable(hex_viewer ${source})
#add_library(hex_viewer ${source})
target_link_libraries(hex_viewer
  zjucad-ptree
  jtf-igl
  jtf-mesh
  AntTweakBar
  glfw3
GLEW
GL
X11
Xxf86vm
Xrandr
Xi
lapack
  )
