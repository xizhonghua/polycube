#! /bin/bash

find -name "*.vtk" > run.log
cat run.log | while read oneline
do
    make -f Makefile VTK_TET_FILE=${oneline}
done