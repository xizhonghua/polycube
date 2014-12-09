#! /bin/bash
HEX_UI=/home/tfjiang/Projects/hex/build/bin/hex_ui
HEX_UTIL=/home/tfjiang/Projects/hex/build/bin/hex_utils
POLYCUBE=/home/tfjiang/Projects/hex/build/bin/polycube

model=$1
size=$2
loop=$3

COUNTER=0
while [ $COUNTER -lt $loop ]; do
echo "${COUNTER}"
#echo The counter is $COUNTER
$POLYCUBE prog=polycube package=hj alg=More output=$model-$size-polycube-split-${COUNTER}.tet linear_solver/type=direct linear_solver/name=cholmod iter=50 tet=$model-$size-split-${COUNTER}.tet
let COUNTER_POST=COUNTER+1
$HEX_UTIL split_not_aligned_face $model-$size-split-${COUNTER}.tet $model-$size-polycube-split-${COUNTER}.tet $model-$size-split-${COUNTER_POST}.tet
let COUNTER=COUNTER+1
done

