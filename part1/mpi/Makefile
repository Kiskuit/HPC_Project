target: decide

CFLAGS=-O3 -std=c99 -g
CC=mpicc

main.o: main.c projet.h
aux.o: aux.c projet.h


decide: main.o aux.o
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean exec

clean:
	rm -f *.o decide ; rm -f *~

exec:
	/usr/bin/time -p mpirun -n 10 -hostfile hostfile --map-by node ./decide "7K//k1P/7p b"
