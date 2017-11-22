CC=/usr/local/gcc-7.1/bin/gcc -O3 -ffast-math
all: fit

fit: fit.c matrixlib.o
	$(CC) -fcilkplus -Wall -std=gnu99 -o fit \
		fit.c matrixlib.o -lm -lcilkrts

matrixlib.o: matrixlib.c matrixlib.h
	$(CC) -fcilkplus -std=gnu99 -c matrixlib.c
