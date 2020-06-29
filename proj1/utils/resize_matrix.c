#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s <input_filename> <resized_rows> <resized_cols> <output_filename>\n", argv[0]);
        return -1;
    }

    int resized_rows = atoi(argv[2]);
    int resized_cols = atoi(argv[3]);
    int rows;
    int cols;
    FILE *infileptr;
    FILE *outfileptr;
    
    infileptr = fopen(argv[1], "r");
    fread(&rows, sizeof(int), 1, infileptr);
    fread(&cols, sizeof(int), 1, infileptr);

    if (rows * cols != resized_rows * resized_cols) {
        printf("Incompatible resize for %d x %d to %d x %d!\n", rows, cols, resized_rows, resized_cols);
        return -1;
    }

    float *storage;
    storage = malloc(rows * cols * sizeof(float));
    fread(storage, sizeof(float), rows * cols, infileptr);

    fclose(infileptr);

    outfileptr = fopen(argv[4], "w");
    fwrite(&resized_rows, sizeof(int), 1, outfileptr);
    fwrite(&resized_cols, sizeof(int), 1, outfileptr);
    fwrite(storage, sizeof(float), rows * cols, outfileptr);

    fclose(outfileptr);
    free(storage);

    return 0;
}