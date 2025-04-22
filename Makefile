#compilers
CC=icpc

LIBS += -L . -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L. -lDC 

all: MAIN

MAIN: MAIN
	$(CC) -I -std=c++11 -qopenmp -O3 main.cpp -Wl,-rpath . $(LIBS) -o main

clean:
	rm -f main
