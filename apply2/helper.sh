#!/bin/bash

mkdir SUPRA_Output
for i in $(seq 0 96)
do
	coord="opt_dir$i/coord"
	opt_struc="opt_dir$i/opt_struc.xyz"
    echo $coord
    echo $opt_struc
	t2x $coord > $opt_struc 2>/dev/null
	mv $opt_struc SUPRA_Output/conformer$i.xyz
	#rm -r opt_dir$i 
done
