#! /bin/bash

#A=refined_80.ext
A=jac_refined_176_0.00_0_ref
for (( i = 60; i < 120; i++)); 
do 
    let "b=i+1"
    echo ${b}
    ~/Projects/hex/build/bin/hex_utils improve_hex ${A}_${i}.hex ${A}_${b}.hex 
    ~/Projects/hex/build/bin/hex_utils hex_jac ${A}_${b}.hex 
done;