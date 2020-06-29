#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>

int samplesInsideCircle(const int numSamples, int seed) {

    int counter = 0;
    double x, y;

    for (int s = 0; s < numSamples; ++s) {
        x = (double) rand_r(&seed) / RAND_MAX;
        y = (double) rand_r(&seed) / RAND_MAX;

        if (x * x + y * y <= 1.0)
            counter++;
    }

    return counter;
}

int main(int argc, char** argv)
{
    int n;
    float pi;
    int seed;
    int threads;
    int global_counter;

    MPI_Init(&argc, &argv);

    int processes;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank) {
        if (argc != 4) {
            printf("Correct command line: ");
            printf("%s <# samples> <# threads> <seed>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    n = atoi(argv[1]);
    threads = atoi(argv[2]);
    seed = atoi(argv[3]);

    omp_set_num_threads(threads);
    int chunk = n / (threads * processes);
    
    int counter = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    double start_time = MPI_Wtime();
    #pragma omp parallel for shared(chunk, seed) reduction(+:counter)
    for (int i = 0; i < threads; ++i)
    {
        counter = samplesInsideCircle(chunk, seed);
    }

    MPI_Reduce(&counter, &global_counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (!rank) {
        pi = 4.0 * global_counter / n;
        printf("Samples: %d, Estimate of pi: %f, Processes: %d, Threads: %d, Using %lf seconds\n", 
                n, pi, processes, threads, (end_time - start_time));
    }

    MPI_Finalize();
    
    return 0;
}
