#!/bin/bash

for folder in opt_dir*
do
	cd $folder
	mkdir gfnff_singlepoint
	cd gfnff_singlepoint
	cp ../struc.xyz .
	xtb --gfnff struc.xyz > gfnff.out
	cd ..
	cd ..
done
