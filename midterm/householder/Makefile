#CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math
CC=gcc
all: fit

fit: fit.c matrixlib.o
	$(CC)  -Wall -std=gnu99 -o fit \
		fit.c matrixlib.o -lm 

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -std=gnu99 -c matrixlib.c
