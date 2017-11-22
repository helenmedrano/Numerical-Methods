#CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math
CC=gcc
all: gramschmidt

gramschmidt: gramschmidt.c matrixlib.o tictoc.o
	$(CC) -Wall -std=gnu99 -o gramschmidt \
		gramschmidt.c matrixlib.o tictoc.o -lm 

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -std=gnu99 -c matrixlib.c

tictoc.o: tictoc.c tictoc.h
	$(CC) -std=gnu99 -c tictoc.c
