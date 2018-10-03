#!/bin/bash
fpack: rnpack.f
	# make clean
	# gfortran -g -fcheck=all -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
	gfortran -g -fcheck=all -std=legacy -o rnpack rnpack.f
	./rnpack

cpack: rnpack.c
	f2c rnpack.f
	g++ -c -w rnpack.c
	g++ -o rnpack rnpack.o -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm -lgfortran -llapacke
	# gcc -o rnpack.o -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm

	# gcc -o rnpack.c -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm
	# gcc -lm -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm rnpack.c 
# *.o rnpack.o



clean:
	rm -rf rnpack
