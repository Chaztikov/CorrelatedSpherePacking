#!/bin/bash

rnpack: rnpack.c
	gcc -o rnpack.o -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm
	# gcc -o rnpack.c -L/usr/lib -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm
	# gcc -lm -I/usr/lib/x86_64-linux-gnu/libf2c.a -lm rnpack.c 
# *.o rnpack.o

# rnpack: rnpack.f
# 	make clean
# 	#gfortran -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
# 	# mpifort -g -fcheck=all -Wall -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
# 	f77 -g -fcheck=all -Wall -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
# 	# mpif90 -g -fcheck=all -Wall -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
# 	# mpif77 -g -fcheck=all -Wall -ffixed-line-length-0 -std=legacy -o rnpack rnpack.f
# 	#gfortran -O -o rnpack rnpack.f
# 	cp rnpack ../EXEC/.
# 	./rnpack

# clean:
# 	rm -rf rnpack
