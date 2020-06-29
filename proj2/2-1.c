#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

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

    if (argc != 4) {
        printf("Correct command line: ");
        printf("%s <# samples> <# threads> <seed>\n", argv[0]);
        return -1;
    }

    n = atoi(argv[1]);
    threads = atoi(argv[2]);
    seed = atoi(argv[3]);

    omp_set_num_threads(threads);
    int chunk = n / threads;
    
    int counter = 0;

    double start_time = omp_get_wtime();
    #pragma omp parallel for shared(chunk, seed) reduction(+:counter)
    for (int i = 0; i < threads; ++i)
    {
        counter = samplesInsideCircle(chunk, seed);
    }
    pi = 4.0 * counter / n;
    double end_time = omp_get_wtime();

    printf("Samples: %d, Estimate of pi: %f, Threads: %d, Using %lf seconds\n", n, pi, threads, (end_time - start_time));
    
    return 0;
}
