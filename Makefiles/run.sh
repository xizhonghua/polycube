#! /bin/bash

rm log
find . -type l > log
MAKEFILE=/home/tfjiang/Projects/hex/2013-siggraph/dense/Makefile_polycube_diff
cat log | while read oneline
do 
	cd ${oneline}
	echo ${oneline}
	mkdir diff-1 diff-2 diff-5
	cd diff-1
	echo "diff-1"
	time make -f ${MAKEFILE} NORMAL_DIFF_W=1 &> log
	cd ../diff-2
	echo "diff-2"
	time make -f ${MAKEFILE} NORMAL_DIFF_W=2 &> log
	cd ../diff-5
	echo "diff-5"
	time make -f ${MAKEFILE} NORMAL_DIFF_W=5 &> log
	cd ../../
done