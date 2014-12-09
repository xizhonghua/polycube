aux_source_directory(. source)
file(GLOB headers *.h)

enable_language (Fortran)
set(SRC_DIR "{${CMAKE_CURRENT_SOURCE_DIR},$ENV{HOME}/.maxima}")

add_custom_command(OUTPUT permutation.f
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b permutation.mac
COMMAND "pwd"
COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i permutation.f
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i permutation.f
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i permutation.f
DEPENDS permutation.mac)

#add_custom_command(OUTPUT yufei_permu.f
#COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b yufei_permu.mac
#COMMAND "pwd"
#COMMAND "sed" ARGS -r "'s/ sqrt\\(([0-9]+)\\)/ sqrt\(\\1\\.0\)/g'" -i yufei_permu.f
#COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i yufei_permu.f
#COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i yufei_permu.f
#DEPENDS yufei_permu.mac)

add_custom_command(OUTPUT Ax.f90
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b Ax.mac
#in maxima 5.30m, still need following two commands
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i Ax.f90
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i Ax.f90
DEPENDS Ax.mac)

add_custom_command(OUTPUT relation.f90
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b relation.mac
#in maxima 5.30m, still need following two commands
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i relation.f90
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i relation.f90
DEPENDS relation.mac)

add_custom_command(OUTPUT rot.f90
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b rot.mac
#in maxima 5.30m, still need following two commands
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i rot.f90
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i rot.f90
DEPENDS rot.mac)

add_custom_command(OUTPUT smooth_L1.f90
COMMAND "maxima" ARGS --userdir=${SRC_DIR} -b smooth_L1.mac
#in maxima 5.30m, still need following two commands
COMMAND "sed" ARGS -r "'s/\\[/\(/g'" -i smooth_L1.f90
COMMAND "sed" ARGS -r "'s/\\]/\)/g'" -i smooth_L1.f90
DEPENDS smooth_L1.mac)


link_directories(
$ENV{HOME}/usr/lib
$ENV{HOME}/usr/lib/Ipopt
)

add_library(hex_frame_opt SHARED ${source} ${headers} permutation.f Ax.f90 relation.f90 rot.f90 smooth_L1.f90)

target_link_libraries(hex_frame_opt
  rot_cubic_f_SH
zjucad-linear-solver
  hj-sparse-util
  hj-function
  common
  tetmesh
  lapack
  jtf-function
   jtf-optimizer
ipopt
coinhsl
coinmetis
  hex_param
)
