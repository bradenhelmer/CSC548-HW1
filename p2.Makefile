SRC := p2.c p2_func.c

all: p2 

p2: $(SRC)
	mpicc -lm -O3 -o p2 p2.c p2_func.c