#!/bin/bash

f="0"
while [ $f -lt 3500 ]; do
	let "n = 250000 + f * 100"
	echo $n
	cd $n/
	rm average_energies_$n
	rm asymmetry_$n
	rm almo_new
	rm energies
	rm almo_energy	
	cd ..
	let "f = f + 1"
done


