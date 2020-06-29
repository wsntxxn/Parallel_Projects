#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void init_matrix(int row, int col, float ***matrix, float **matrix_storage) {
    // float **matrix = (float **) malloc(sizeof(float *) * row);
    *matrix = (float **) malloc(sizeof(float *) * row);
    *matrix_storage = (float *) malloc(sizeof(float) * row * col);
    float *ptr = *matrix_storage;
    for (int i = 0; i < row; ++i) {
        // matrix[i] = (float*) malloc(sizeof(float) * col);
        *(*matrix + i) = ptr;
        ptr += col;
    }
    
}

void im2row(float **data_img, 
            int width, 
            int height,
            int ksize,
            int stride, 
            int *width_out, 
            int *height_out,
            float ***data_row,
            float **data_row_storage) {

    *width_out = (width - ksize) / stride + 1;
    *height_out = (height - ksize) / stride + 1;
    // printf("out_w: %d out_h: %d\n", width_out, height_out);
    init_matrix(*width_out * *height_out, ksize * ksize, data_row, data_row_storage);

    for (int i = 0; i < *height_out; ++i) {

        for (int j = 0; j < *width_out; ++j) {
            // assign [i][j] position of new data

            for (int h = 0; h < ksize; ++h) {
                for (int w = 0; w < ksize; ++w) {
                    // printf("data_img: %f\n", data_img[i + h][j + w]);
                    (*data_row)[i * (*width_out) + j][h * ksize + w] = data_img[i * stride + h][j * stride + w];
                }
            }
        }
    }
}

void write_matrix(int row, int col, float **matrix_storage, char* filename){
    FILE *fout = fopen(filename, "w");
    fwrite(&row, sizeof(int), 1, fout);
    fwrite(&col, sizeof(int), 1, fout);
    fwrite(*matrix_storage, sizeof(float), row * col, fout);

    fclose(fout);
}

int main(int argc, char* argv[]){
    float **matrix;
    float *matrix_storage;

    if (argc < 7) {
        printf("Usage: %s <rows> <cols> <ksize> <stride> <filename> <is kernel> [<raw matrix filename> <is pooling>]\n", argv[0]);
        return -1;
    }

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int ksize = atoi(argv[3]);
    int stride = atoi(argv[4]);
    int is_kernel = atoi(argv[6]);
    int is_pooling;

    if (is_kernel)
        is_pooling = atoi(argv[8]);

    srand((unsigned) time(NULL));

    init_matrix(m, n, &matrix, &matrix_storage);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (is_pooling)
                matrix[i][j] = 1.0 / (ksize * ksize);
            else {
                int sign = ((float) rand() / RAND_MAX > 0.5? 1: -1); 
                matrix[i][j] = (float) 2 * rand() / RAND_MAX * sign;
            }
        }
    }
    
    /*
    printf("initial matrix: \n");
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            printf("%.2f\t", matrix[i][j]);
        printf("\n");
    }
    printf("\n");
    */

    if (is_kernel) {
        write_matrix(m, n, &matrix_storage, argv[5]);
        free(matrix_storage);
        free(matrix);
        return 0;
    }
    else {
        write_matrix(m, n, &matrix_storage, argv[7]);
    }

    int width_out;
    int height_out;
    float **matrix_row;
    float *matrix_row_storage;
    im2row(matrix, n, m, ksize, stride, &width_out, &height_out, &matrix_row, &matrix_row_storage);
    printf("Output row: %d, Output column: %d\n", width_out, height_out);

    int row_out = height_out * width_out;
    int col_out = ksize * ksize;

    /*
    printf("im2row transformed matrix: \n");
    for (int i = 0; i < row_out; ++i) {
        for (int j = 0; j < col_out; ++j)
            printf("%.2f\t", matrix_row[i][j]);
        printf("\n");
    }
    */

    FILE *fout = fopen(argv[5], "w");
    fwrite(&row_out, sizeof(int), 1, fout);
    fwrite(&col_out, sizeof(int), 1, fout);
    fwrite(matrix_row_storage, sizeof(float), row_out * col_out, fout);

    fclose(fout);
    free(matrix_storage);
    free(matrix);
    free(matrix_row_storage);
    free(matrix_row);
    return 0;
}