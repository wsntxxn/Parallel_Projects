MPICC?=mpicc
MPICXX?=mpicxx
ALL=1-1 1-2 1-3-small 1-3-big

all: ${ALL}

1-1: 1-1.c
	${MPICC} -o $@ 1-1.c

1-2: 1-2.c
	${MPICC} -o $@ 1-2.c -lm

1-3-small: 1-3-small.cpp
	${MPICXX} -o $@ 1-3-small.cpp

1-3-big: 1-3-big.cpp
	${MPICXX} -o $@ 1-3-big.cpp

clean:
	rm -f ${ALL}
