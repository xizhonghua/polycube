file(GLOB header *.h)
file(GLOB source *.cpp)

enable_language (Fortran)

set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/.maxima}")

add_custom_command(OUTPUT tri-normal.f
   COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b tri-normal.mac
   COMMAND "pwd"
   COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i tri-normal.f
   COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i tri-normal.f
   COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i tri-normal.f
   DEPENDS tri-normal.mac)

add_custom_command(OUTPUT tri-area.f
   COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b tri-area.mac
   COMMAND "pwd"
   COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i tri-area.f
   COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i tri-area.f
   COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i tri-area.f
   DEPENDS tri-area.mac)

add_custom_command(OUTPUT tri-area-normal.f
   COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b tri-area-normal.mac
   COMMAND "pwd"
   COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i tri-area-normal.f
   COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i tri-area-normal.f
   COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i tri-area-normal.f
   DEPENDS tri-area-normal.mac)

add_custom_command(OUTPUT so3.f
  COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b so3.mac
  COMMAND "pwd"
  COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i so3.f
  COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i so3.f
  COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i so3.f
  DEPENDS so3.mac)

add_custom_command(OUTPUT tri_area_normal_diff.f90
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b tri_area_normal_diff.mac
#in maxima 5.30m, still need following two commands
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i tri_area_normal_diff.f90
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i tri_area_normal_diff.f90
DEPENDS tri_area_normal_diff.mac)

add_custom_command(OUTPUT tet-vol.f90
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b tet-vol.mac
#in maxima 5.30m, still need following two commands
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i tet-vol.f90
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i tet-vol.f90
DEPENDS tet-vol.mac)

add_library(mesh_func SHARED ${header} ${source} tri-normal.f tri-area.f tri-area-normal.f tet-vol.f90
so3.f tri_area_normal_diff.f90)
