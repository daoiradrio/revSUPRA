#!/bin/bash

mkdir workdir

for struc in conformer*
do
	mv $struc workdir
	x2t workdir/$struc > workdir/coord
	cp control workdir
	cd workdir
	uff > uff.out 2>&1
	cd ..
	t2x workdir/coord > workdir/$struc
	mv workdir/$struc $struc
	rm -f workdir/*	
done

rm -r workdir
