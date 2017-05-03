target: decide

CFLAGS=-O3 -std=c99 -g -fdiagnostics-color=auto
CC=mpicc

main.o: main.c projet.h
aux.o: aux.c projet.h


decide: main.o aux.o
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f *.o decide ; rm -f *~

