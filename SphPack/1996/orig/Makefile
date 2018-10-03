#!/bin/bash
rnpack: rnpack.f
	make clean
	#gfortran -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
	gfortran -g -fcheck=all -Wall -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
	#gfortran -O -o rnpack rnpack.f
	cp rnpack ../EXEC/.
	./rnpack
clean:
	rm -rf rnpack
