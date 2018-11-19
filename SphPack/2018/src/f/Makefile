#!/bin/bash
fpack: rnpack.f
	# make clean
	# gfortran -g -fcheck=all -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
	gfortran -g -lgfortran -fcheck=all -std=legacy -o rnpack rnpack.f
	./rnpack

cpack: rnpack.c
	f2c rnpack.f
	icc -c -w rnpack.c
	icc -o rnpack rnpack.o main_testrnpack.c -L/usr/lib -I/usr/include 
	#-L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a
	# gcc -o rnpack.o -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm

	# gcc -o rnpack.c -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm
	# gcc -lm -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm rnpack.c 
# *.o rnpack.o



clean:
	rm -rf rnpack
