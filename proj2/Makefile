all: 2-1 2-1-MPI 2-2 2-2-MPI 2-3 2-3-MPI

2-1: 2-1.c
	gcc --openmp -o 2-1 2-1.c

2-1-MPI: 2-1-MPI.c
	mpicc --openmp -o 2-1-MPI 2-1-MPI.c

2-2: 2-2.c
	gcc --openmp -o 2-2 2-2.c

2-2-MPI: 2-2-MPI.c
	mpicc --openmp -o 2-2-MPI 2-2-MPI.c

2-3: 2-3.c
	gcc --openmp -o 2-3 2-3.c

2-3-MPI: 2-3-MPI.c
	mpicc --openmp -o 2-3-MPI 2-3-MPI.c

clean:
	rm 2-1 2-1-MPI 2-2 2-2-MPI 2-3 2-3-MPI
