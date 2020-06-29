SEED=1234

echo -e "Running 2-1 (openmp)"
./2-1 10000000 4 $SEED
echo -e "Running 2-1 (openmp + MPI)"
mpirun -np 2 -f ./config ./2-1-MPI 10000000 2 $SEED

echo -e "Running 2-2 (openmp)"
./2-2 4
echo -e "Running 2-2 (openmp + MPI)"
mpirun -np 2 -f ./config ./2-2-MPI 2

echo -e "Running 2-3 (openmp)"
./2-3 4 $SEED
echo -e "Running 2-3 (openmp + MPI)"
mpirun -np 2 -f ./config ./2-3-MPI 2 $SEED
