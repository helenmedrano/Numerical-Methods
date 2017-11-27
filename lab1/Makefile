all: invpower

clean:
	rm -f p2q1 p2q1.eps p2q2 ec *.o

invpower: p2q2.c matrixlib.o
	gcc -Wall --std=gnu99 -o invpower p2q2.c matrixlib.o -lm

matrixlib.o: matrixlib.c matrixlib.h
	gcc -Wall --std=gnu99 -c matrixlib.c
