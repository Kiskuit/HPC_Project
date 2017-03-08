target: decide

CFLAGS=-g -O3 -std=c99 -fopenmp -fdiagnostics-color=auto

main.o: main.c projet.h
aux.o: aux.c projet.h


decide: main.o aux.o
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f *.o decide

