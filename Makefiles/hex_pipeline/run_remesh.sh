#! /bin/bash

TIME_NOW=`date --iso=second`

HexUtil=/home/$(whoami)/Projects/hex/build/bin/hex_utils

mkdir fix
cd fix
time make -f ../Makefile_fix_singularity &> log_${TIME_NOW}
cd ../

mkdir equation-no-surface
cd equation-no-surface
time make -f ../Makefile_equation-no-surface &> log_${TIME_NOW}
cd ../

mkdir param-no-surface
cd param-no-surface
#time make -f ../Makefile_remesh2 &> log_${TIME_NOW}
time make -f ../Makefile_remesh &> log_${TIME_NOW}
cd ..

mkdir equation-with-l1-surface-type
cd equation-with-l1-surface-type
time make -f ../Makefile_equation-with-surface &> log_${TIME_NOW}
cd ..

mkdir param-with-l1-surface-type
cd param-with-l1-surface-type
time make -f ../Makefile_fix_remesh &> log_${TIME_NOW}
cd ..

mkdir hex
cd hex
${HexUtil} tet2vtk ../param-with-l1-surface-type/output.integer.tet &> param.tet.vtk
${HexUtil} tet2vtk ../param-no-surface/after_l1.uncut_tet.tet &> orig.tet.vtk 
ln -s ../param-no-surface/after_l1.inner_type inner_type
time make -f ../Makefile_hex &> log_${TIME_NOW}