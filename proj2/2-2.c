#include <assert.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

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


int main(int argc, char *argv[]) {	
    int threads;
	if (argc < 2) {
        printf("Correct command line: ");
		printf("%s <# threads> \n", argv[0]);
        return -1;
	}
	int input_size = 1e7; // Last must be 0
    threads = atoi(argv[1]);

    omp_set_num_threads(threads);

	// Set seed to get always same numbers
	srandom(1);

    int *data = generate_data(input_size);

    double start_time = omp_get_wtime();
    #pragma omp parallel
    {
        #pragma omp single
        {
            quickSort(data, 0, input_size - 1);
        }
    }
    double end_time = omp_get_wtime();
    verify(data, input_size);
    free(data);

    printf("Threads: %d, Using %f secs\n", threads, end_time - start_time);
	
	return 0;
}

