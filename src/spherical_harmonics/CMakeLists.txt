enable_language (Fortran)

set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/.maxima}")
message(${SRC_DIR})

add_custom_command(OUTPUT rot_cubic_f_SH.lisp
  COMMAND "maxima" ARGS --userdir="${SRC_DIR}" -b rot_cubic_f_SH.mac
  DEPENDS rot_cubic_f_SH.mac SH.mac)

add_custom_command(OUTPUT rot_cubic_f_SH.f
  COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b rot_cubic_f_SH_fortran.mac
  COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i rot_cubic_f_SH.f
  COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i rot_cubic_f_SH.f
  COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i rot_cubic_f_SH.f
  DEPENDS rot_cubic_f_SH.lisp rot_cubic_f_SH_fortran.mac)

add_library(rot_cubic_f_SH SHARED rot_cubic_f_SH.f sh2zyz.cpp)

target_link_libraries(rot_cubic_f_SH
  minpack)
