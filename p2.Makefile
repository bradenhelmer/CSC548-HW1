SRC := p2.c p2_func.c

all: clean p2

p2: $(SRC)
	mpicc -lm -O3 -o p2 p2.c p2_func.c

serial: p2_serial.c
	mpicc -lm -O3 -o p2_serial p2_serial.c p2_func.c

run: p2
	prun ./p2 100 0 0

clean:
	rm -f p2 *.dat *.png
