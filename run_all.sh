#!/bin/bash

f="0"
while [ $f -lt 3500 ]; do
	let "n = 250000 + f * 100"
	echo $n
	cp ./almo_new $n/
	cd $n/
	rm average_energies_$n
	rm asymmetry_$n
	./almo_new
	mv average_energies average_energies_$n
	mv asymmetry asymmetry_$n
	cp average_energies_$n ../all_energies
	cp asymmetry_$n ../all_asymmetries
	cd ..
	let "f = f + 1"
done


