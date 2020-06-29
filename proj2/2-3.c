#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define MAX_OUT_DEGREE 10
#define GRAPH_SIZE 1024000
#define ITERATION 100
#define PRINT_EVERY 10

int duplicate (int index, int* indexes, int size) {
    for (int i = 0; i < size; ++i) {
        if (indexes[i] == index)
            return 1;
    }
    return 0;
}

typedef struct Node {
    int index;
    int out_degree;
    int in_degree;
    int *parents;
} Node;

typedef struct Graph {
    int size;
    struct Node *nodes;
} Graph;

struct Graph* create_graph (int size) {
    Graph *graph;
    graph = (Graph*) malloc(sizeof(Graph));
    graph -> size = size;
    graph -> nodes = (Node *) malloc(sizeof(Node) * size);

    for (int i = 0; i < size; ++i) {
        graph -> nodes[i].index = i;
        graph -> nodes[i].in_degree = 0;
        graph -> nodes[i].parents = NULL;
    }

    // create graph
    for (int i = 0; i < size; ++i) {
        int out_degree = rand() % MAX_OUT_DEGREE + 1;
        graph -> nodes[i].out_degree = out_degree;
        
        int out_node_indexes[MAX_OUT_DEGREE];
        for (int j = 0; j < out_degree; ++j) {
            int out_node_index = rand() % size;
            while (duplicate(out_node_index, out_node_indexes, j)) {
                out_node_index = rand() % size;
            }
            out_node_indexes[j] = out_node_index;

            graph -> nodes[out_node_index].in_degree++;
            int in_degree = graph -> nodes[out_node_index].in_degree;
            graph -> nodes[out_node_index].parents = realloc(graph -> nodes[out_node_index].parents, sizeof(int) * in_degree);
            graph -> nodes[out_node_index].parents[in_degree - 1] = i;
        }
        /*printf("%d %d\n", i, graph -> nodes[17].parents[0].index);*/
    }

    return graph;
}

int main(int argc, char* argv[]) {
    
    if (argc != 3) {
        printf("Correct command line: ");
        printf("%s <# threads> <seed>\n", argv[0]);
        return -1;
    }

    int threads = atoi(argv[1]);
    omp_set_num_threads(threads);
    int seed = atoi(argv[2]);
    srand(seed);
    Graph* graph;
    graph = create_graph(GRAPH_SIZE);

    float *v;
    float *v_old;
    
    v = (float *) malloc(sizeof(float) * GRAPH_SIZE);
    v_old = (float *) malloc(sizeof(float) * GRAPH_SIZE);

    double start_time = omp_get_wtime();

    for (int i = 0; i < GRAPH_SIZE; ++i)
        v_old[i] = 1.0 / GRAPH_SIZE;
    

    // pagerank 
    for (int i = 0; i < ITERATION; ++i) {

        #pragma omp parallel for
        // calculate new value
        for (int j = 0; j < GRAPH_SIZE; ++j) {
            v[j] = 0.0;
            for (int k = 0; k < graph -> nodes[j].in_degree; ++k) {
                v[j] += v_old[graph -> nodes[j].parents[k]] / graph -> nodes[graph -> nodes[j].parents[k]].out_degree;
            }
        }

        if (i % PRINT_EVERY == 0) {
            // calculate error
            float error = 0;
            #pragma omp parallel for
            for (int j = 0; j < GRAPH_SIZE; ++j)
                error += v[j] - v_old[j] > 0? v[j] - v_old[j]: v_old[j] - v[j];
                printf("Iter %d, Error: %f\n", i, error); 
        }

        #pragma omp parallel for
        // update v
        for (int j = 0; j < GRAPH_SIZE; ++j)
            v_old[j] = v[j];
    }


    double end_time = omp_get_wtime();
    printf("Threads: %d Using %lf seconds\n", threads, (end_time - start_time));


    for (int i = 0; i < GRAPH_SIZE; ++i) {
        free(graph -> nodes[i].parents);
    }
    free(graph -> nodes);

    return 0;
}
