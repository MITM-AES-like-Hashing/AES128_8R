CC= g++
CFLAGS= -fopenmp -maes -mavx2 -O3 -std=c++11
DFLAGS= -Domp_nb_threads=16 -DDoFF=8

all: F32 F40 F48 F56 F64

F32:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -DPARTIAL=32

F40:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -DPARTIAL=40

F48:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -DPARTIAL=48

F56:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -DPARTIAL=56

F64:
	$(CC) -o $@ main.cpp $(CFLAGS) $(DFLAGS) -DPARTIAL=64

.PHONY:clean
clean:
	rm F32 F40 F48 F56 F64