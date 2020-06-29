#include <assert.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int* generate_data(int n)
{
	int i=0;
	int *data = malloc(sizeof(int) * n);
	for (i = 0; i < n; i++)
		data[i] = random();
	
	return data;
}

int Partition(int *data, int low, int high) {
    int pivot = data[low];
    while (low < high) {
        while (low < high && data[high] >= pivot) high--;
        data[low] = data[high];
        while (low < high && data[low] <= pivot) low++;
        data[high] = data[low];
    }
    data[low] = pivot;
    return low;
}

void quickSort(int *data, int low, int high)
{
    
    if (low < high) {
        int split = Partition(data, low, high);
        #pragma omp task
        quickSort(data, low, split - 1);
        #pragma omp task
        quickSort(data, split + 1, high);
    }
}


void verify(int *data, int n)
{
	int i;
	for (i = 1; i < n; i++)
		assert(data[i] >= data[i-1]);
}

void merge_array(int *data, int length, int chunks) {
    int chunk_length = length / chunks;
    int *cur_chunk_idxs = malloc(sizeof(int) * chunks);
    for (int i = 0; i < chunks; ++i)
        cur_chunk_idxs[i] = i * chunk_length;

    int *data_tmp = malloc(sizeof(int) * length);

    int cur_idx = 0;

    while (cur_idx < length) {
        int cur_min = 0x7fffffff;
        int cur_min_chunk_idx;
        for (int j = 0; j < chunks; ++j) {
            if (cur_chunk_idxs[j] < (j + 1) * chunk_length &&
                    data[cur_chunk_idxs[j]] < cur_min) {
                cur_min = data[cur_chunk_idxs[j]];
                cur_min_chunk_idx = j;
            }
        }
        data_tmp[cur_idx] = cur_min;
        cur_idx++;
        cur_chunk_idxs[cur_min_chunk_idx]++;
    }
    
    for (int i = 0; i < length; ++i)
        data[i] = data_tmp[i];
    
    free(cur_chunk_idxs);
    free(data_tmp);
}

int main(int argc, char *argv[]) {	
    int threads;
    int processes;
    int rank;
    int *data;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!rank) {
        if (argc < 2) {
            printf("Correct command line: ");
            printf("%s <# threads> \n", argv[0]);
            return -1;
        }
    }

	int input_size = 1e7;
    threads = atoi(argv[1]);

    omp_set_num_threads(threads);

	// Set seed to get always same numbers
	srandom(1);

    if (!rank)
        data = generate_data(input_size);

    int send_count = input_size / processes;
    int *local_data;
    local_data = malloc(sizeof(int) * send_count);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    MPI_Scatter(data, send_count, MPI_INT, local_data, send_count, MPI_INT, 0, MPI_COMM_WORLD);

    #pragma omp parallel
    {
        #pragma omp single
        {
            quickSort(local_data, 0, send_count - 1);
        }
    }
    
    MPI_Gather(local_data, send_count, MPI_INT, data, send_count, MPI_INT, 0, MPI_COMM_WORLD);
    free(local_data);
    
    if (!rank) {
        merge_array(data, input_size, processes);
        double end_time = MPI_Wtime();
        verify(data, input_size);
        printf("Processes: %d, Threads: %d, Using %f secs\n", 
                processes, threads, end_time - start_time);
        free(data);
    }
    
    MPI_Finalize();
	
	return 0;
}

