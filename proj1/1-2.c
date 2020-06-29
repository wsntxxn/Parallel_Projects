/*
 * 
 * Matrix multiplication using Cannon's Algorithm in MPI
 *
 * The program takes three command-line arguments: fileA, fileB, and
 * fileC. The first two files contain matrix A and B as the input. The
 * third file is used to store the result matrix C as the output. The
 * program compute: C = A x B. The program assumes the number of 
 * processors p is square.
 *
 * The files containing the matrices are all binary files and have the
 * following format. The matrix is stored in row-wise order and
 * preceded with two integers that specify the dimensions of the
 * matrix. The matrix elements are floating point numbers.
 *
 */

#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

/* set this parameter to reflect the cache line size of the particular
   machine you're running this program */
#define CACHE_SIZE 32768

/* in case later we decide to use another data type */
#define mpitype MPI_FLOAT
typedef float datatype;

/* block decomposition macros */
#define BLOCK_LOW(id, p, n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j, p, n) (((p)*((j)+1)-1)/(n))

int cfileexists(const char *filename) {
    /* try to open file to read */
    FILE *file;
    if (file = fopen(filename, "r")) {
        fclose(file);
        return 1;
    }
    return 0;
}

/* print out error message and exit the program */
void my_abort(const char *fmt, ...) {
    int id;     /* process rank */
    va_list ap; /* argument list */

    va_start(ap, fmt);

    /* only process 0 reports */
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (!id) vprintf(fmt, ap);

    va_end(ap);

    /* all MPI processes exit at this point */
    exit(1);
}

void sum(int rows, int cols, datatype *const *c, datatype *const *partial_c_matrix) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            partial_c_matrix[i][j] += c[i][j];
        }
    }
}

void check_input_files(char *const *argv) {
    if (cfileexists(argv[1]) == 0) {
        my_abort("Error: File %s is not accessible\n", argv[1]);
    }

    if (cfileexists(argv[2]) == 0) {
        my_abort("Error: File %s is not accessible\n", argv[2]);
    }
}

datatype **init_partial_c_matrix(int rows, int cols) {
    datatype **partial_c_matrix = calloc((size_t) rows, sizeof(datatype *));

    for (int i = 0; i < rows; ++i) {
        partial_c_matrix[i] = calloc((size_t) cols, sizeof(datatype));

        for (int j = 0; j < cols; ++j) {
            partial_c_matrix[i][j] = 0;
        }
    }
    return partial_c_matrix;
}

datatype **init_local_c(int rows, int cols, datatype **c, datatype *sc) {
    sc = (datatype *) malloc(rows * cols * sizeof(datatype));
    memset(sc, 0, rows * cols * sizeof(datatype));

    c = (datatype **) malloc(rows * sizeof(datatype *));

    for (int j = 0; j < rows; j++) {
        c[j] = &sc[j * cols];
    }

    return c;
}


/* return the data size in bytes */
int get_size(MPI_Datatype t) {
    if (t == MPI_BYTE) return sizeof(char);
    else if (t == MPI_DOUBLE) return sizeof(double);
    else if (t == MPI_FLOAT) return sizeof(float);
    else if (t == MPI_INT) return sizeof(int);
    else {
        printf("Error: Unrecognized argument to 'get_size'\n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -3);
    }
    return 0;
}


/* allocate memory from heap */
void *my_malloc(int id, int bytes) {
    void *buffer;
    if ((buffer = malloc((size_t) bytes)) == NULL) {
        printf("Error: Malloc failed for process %d\n", id);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -2);
    }
    return buffer;
}

/* Read a matrix from a file. */
void read_checkerboard_matrix(
        char *s,              /* IN - File name */
        void ***subs,         /* OUT - 2D array */
        void **storage,       /* OUT - Array elements */
        MPI_Datatype dtype,   /* IN - Element type */
        int *rows,            /* OUT - Array rows */
        int *cols,            /* OUT - Array cols */
        MPI_Comm grid_comm,   /* IN - Communicator */
        int max_dim           /* IN - Dimenson holding maximum (last) elements*/)
{
    void *buffer;         /* File buffer */
    int coords[2];      /* Coords of proc receiving
                                 next row of matrix */
    int datum_size;     /* Bytes per elements */
    int dest_id;        /* Rank of receiving proc */
    int grid_coord[2];  /* Process coords */
    int grid_id;        /* Process rank */
    int grid_period[2]; /* Wraparound */
    int grid_size[2];   /* Dimensions of grid */
    int i, j, k;
    FILE *infileptr;      /* Input file pointer */
    void *laddr;          /* Used when proc 0 gets row */
    int local_cols;     /* Matrix cols on this proc */
    int local_rows;     /* Matrix rows on this proc */
    int max_rows;       /* Maximum matrix rows in this communicator */
    int max_cols;       /* Maximum matrix cols in this communicator */
    void **lptr;           /* Pointer into 'subs' */
    int p;              /* Number of processes */
    void *raddr;          /* Address of first element
                                 to send */
    void *rptr;           /* Pointer into 'storage' */
    MPI_Status status;         /* Results of read */

    MPI_Comm_rank(grid_comm, &grid_id);
    MPI_Comm_size(grid_comm, &p);
    datum_size = get_size(dtype);

    /* Process 0 opens file, gets number of rows and
       number of cols, and broadcasts this information
       to the other processes. */

    if (grid_id == 0) {
        infileptr = fopen(s, "r");
        if (infileptr == NULL) *rows = 0;
        else {
            fread(rows, sizeof(int), 1, infileptr);
            fread(cols, sizeof(int), 1, infileptr);
        }
    }
    MPI_Bcast(rows, 1, MPI_INT, 0, grid_comm);

    if (!(*rows)) MPI_Abort(MPI_COMM_WORLD, -1);

    MPI_Bcast(cols, 1, MPI_INT, 0, grid_comm);

    /* Each process determines the size of the submatrix
       it is responsible for. */

    MPI_Cart_get(grid_comm, 2, grid_size, grid_period,
                 grid_coord);
    local_rows = BLOCK_SIZE(grid_coord[0], grid_size[0], *rows);
    local_cols = BLOCK_SIZE(grid_coord[1], grid_size[1], *cols);
    max_rows = BLOCK_SIZE(grid_size[0] - 1, grid_size[0], *rows);
    max_cols = BLOCK_SIZE(grid_size[1] - 1, grid_size[1], *cols);


    /* Dynamically allocate two-dimensional matrix 'subs' */

    if (max_dim == 0) { /* for B, holding maximum row */
        *storage = my_malloc(grid_id,
                             max_rows * local_cols * datum_size);
        memset(*storage, 0, max_rows * local_cols * datum_size);
        *subs = (void **) my_malloc(grid_id, max_rows * sizeof(void *));
        lptr = (void *) *subs;
        rptr = (void *) *storage;
        for (i = 0; i < max_rows; i++) {
            *(lptr++) = (void *) rptr;
            rptr += local_cols * datum_size;
        }

    }
    else {
        *storage = my_malloc(grid_id,
                             local_rows * max_cols * datum_size);
        memset(*storage, 0, local_rows * max_cols * datum_size);
        *subs = (void **) my_malloc(grid_id, local_rows * sizeof(void *));
        lptr = (void *) *subs;
        rptr = (void *) *storage;
        for (i = 0; i < local_rows; i++) {
            *(lptr++) = (void *) rptr;
            rptr += max_cols * datum_size;
        }
    }

    /* Grid process 0 reads in the matrix one row at a time
       and distributes each row among the MPI processes. */

    if (grid_id == 0)
        buffer = my_malloc(grid_id, *cols * datum_size);

    /* For each row of processes in the process grid... */
    for (i = 0; i < grid_size[0]; i++) {
        coords[0] = i;

        /* For each matrix row controlled by this proc row...*/
        for (j = 0; j < BLOCK_SIZE(i, grid_size[0], *rows); j++) {

            /* Read in a row of the matrix */

            if (grid_id == 0) {
                fread(buffer, datum_size, *cols, infileptr);
            }

            /* Distribute it among process in the grid row */

            for (k = 0; k < grid_size[1]; k++) {
                coords[1] = k;

                /* Find address of first element to send */
                raddr = buffer +
                        BLOCK_LOW(k, grid_size[1], *cols) * datum_size;

                /* Determine the grid ID of the process getting
                   the subrow */
                MPI_Cart_rank(grid_comm, coords, &dest_id);

                /* Process 0 is responsible for sending...*/
                if (grid_id == 0) {

                    /* It is sending (copying) to itself */
                    if (dest_id == 0) {
                        laddr = (*subs)[j];
                        memcpy (laddr, raddr,
                                local_cols * datum_size);

                        /* It is sending to another process */
                    } else {
                        MPI_Send(raddr,
                                 BLOCK_SIZE(k, grid_size[1], *cols), dtype,
                                 dest_id, 0, grid_comm);
                    }

                    /* Process 'dest_id' is responsible for
                       receiving... */
                } else if (grid_id == dest_id) {
                    MPI_Recv((*subs)[j], local_cols, dtype, 0,
                             0, grid_comm, &status);
                }
            }
        }
    }
    if (grid_id == 0) free(buffer);
}

/*
 * Write a matrix distributed in checkerboard fashion to a file.
 */
void write_checkerboard_matrix(
        char *s,                /* IN -File name */
        void **a,               /* IN -2D matrix */
        MPI_Datatype dtype,     /* IN -Matrix element type */
        int m,                  /* IN -Matrix rows */
        int n,                  /* IN -Matrix columns */
        MPI_Comm grid_comm)     /* IN -Communicator */
{
    void *buffer;         /* Room to hold 1 matrix row */
    int coords[2];      /* Grid coords of process
                                 sending elements */
    int datum_size;     /* Bytes per matrix element */
    int els;            /* Elements received */
    int grid_coords[2]; /* Coords of this process */
    int grid_id;        /* Process rank in grid */
    int grid_period[2]; /* Wraparound */
    int grid_size[2];   /* Dims of process grid */
    int i, j, k;
    void *laddr;          /* Where to put subrow */
    int local_cols;     /* Matrix cols on this proc */
    int p;              /* Number of processes */
    int src;            /* ID of proc with subrow */
    MPI_Status status;         /* Result of receive */
    FILE *outfileptr;     /* Output file */

    MPI_Comm_rank(grid_comm, &grid_id);
    MPI_Comm_size(grid_comm, &p);
    datum_size = get_size(dtype);

    if (grid_id == 0) {
        outfileptr = fopen(s, "w");
        if (!outfileptr ||
            fwrite(&m, sizeof(int), 1, outfileptr) != 1 ||
            fwrite(&n, sizeof(int), 1, outfileptr) != 1)
            MPI_Abort(MPI_COMM_WORLD, -1);
    }

    MPI_Cart_get(grid_comm, 2, grid_size, grid_period,
                 grid_coords);
    local_cols = BLOCK_SIZE(grid_coords[1], grid_size[1], n);

    if (!grid_id)
        buffer = my_malloc(grid_id, n * datum_size);

    /* For each row of the process grid */
    for (i = 0; i < grid_size[0]; i++) {
        coords[0] = i;

        /* For each matrix row controlled by the process row */
        for (j = 0; j < BLOCK_SIZE(i, grid_size[0], m); j++) {

            /* Collect the matrix row on grid process 0 and
               print it */
            if (!grid_id) {
                for (k = 0; k < grid_size[1]; k++) {
                    coords[1] = k;
                    MPI_Cart_rank(grid_comm, coords, &src);
                    els = BLOCK_SIZE(k, grid_size[1], n);
                    laddr = buffer +
                            BLOCK_LOW(k, grid_size[1], n) * datum_size;
                    if (src == 0) {
                        memcpy (laddr, a[j], els * datum_size);
                    } else {
                        MPI_Recv(laddr, els, dtype, src, 0,
                                 grid_comm, &status);
                    }
                }

                if (fwrite(buffer, datum_size, n, outfileptr) != n)
                    MPI_Abort(MPI_COMM_WORLD, -1);
            } else if (grid_coords[0] == i) {
                MPI_Send(a[j], local_cols, dtype, 0, 0,
                         grid_comm);
            }
        }
    }
    if (!grid_id) {
        free(buffer);
        fclose(outfileptr);
    }
}

void matmul(datatype **a, datatype **b, datatype **c, int m, int n, int k) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < k; ++j) {
            c[i][j] = 0;
            for (int l = 0; l < n; ++l) {
                c[i][j] += a[i][l] * b[l][j];
            }
        }
    }

}


int main(int argc, char *argv[]) {
    int p, p_sqrt;
    int id, coord[2];
    int dim[2], period[2];
    MPI_Comm comm;
    int ma, na, mb, nb, local_rows, local_cols, local_a_col_b_row, max_a_col_b_row;
    datatype **a, *sa;
    datatype **b, *sb;
    datatype **c, *sc;

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    /* start couting time */
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed_time = -MPI_Wtime();

    /* make sure the number of arguments is correct */
    if (argc != 4) {
        my_abort("Usage: %s fileA fileB fileC\n", argv[0]);
    }

    /* create 2D cartesion communicator and obtain the system configurations */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    p_sqrt = (int) sqrt(p);

    if (p_sqrt * p_sqrt != p) {
        my_abort("Error: number of processors (p=%d) must be a square number\n", p);
    }

    dim[0] = dim[1] = p_sqrt;
    period[0] = period[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, 0, &comm);
    MPI_Comm_rank(comm, &id);
    MPI_Cart_coords(comm, id, 2, coord);

    /* checking input files are accessible */
    check_input_files(argv);

    /* read the submatrix of A managed by this process */
    read_checkerboard_matrix(argv[1], (void ***) &a, (void **) &sa, mpitype, &ma, &na, comm, 1);

    /*if (id == 2) {*/
        /*int local_rows = BLOCK_SIZE(coord[0], dim[0], ma);*/
        /*int max_cols = BLOCK_SIZE(dim[1] - 1, dim[1], na);*/
        /*for (int i = 0; i < local_rows; ++i) {*/
            /*for (int j = 0; j < max_cols; ++j)*/
                /*printf("%.2lf\t", a[i][j]);*/
            /*printf("\n");*/
        /*}*/
    /*}*/

    /* read the submatrix of B managed by this process */
    read_checkerboard_matrix(argv[2], (void ***) &b, (void **) &sb, mpitype, &mb, &nb, comm, 0);
    /*if (id == 3) {*/
        /*int max_rows = BLOCK_SIZE(dim[0] - 1, dim[0], mb);*/
        /*int local_cols = BLOCK_SIZE(coord[1], dim[1], nb);*/
        /*for (int i = 0; i < max_rows; ++i) {*/
            /*for (int j = 0; j < local_cols; ++j)*/
                /*printf("%.2lf\t", b[i][j]);*/
            /*printf("\n");*/
        /*}*/
    /*}*/

    if (na != mb) {
        my_abort("id = %d, colume of matrix A does not equal row of B\n", id);
    }

    /* IMPORTANT: we don't have the entire matrix; only the sub */
    local_rows = BLOCK_SIZE(coord[0], dim[0], ma);
    local_cols = BLOCK_SIZE(coord[1], dim[1], nb);
    local_a_col_b_row = BLOCK_SIZE(coord[1], dim[1], na);
    max_a_col_b_row = BLOCK_SIZE(dim[1] - 1, dim[1], na);

    int source, dest;

    datatype **partial_c_matrix = init_partial_c_matrix(local_rows, local_cols);

    MPI_Cart_shift(comm, 1, -coord[0], &source, &dest);
    MPI_Sendrecv_replace(sa, local_rows * max_a_col_b_row, mpitype, dest, 0, source, 0, comm, MPI_STATUS_IGNORE);

    MPI_Sendrecv_replace(&local_a_col_b_row, 1, MPI_INT, dest, 0, source, 0, comm, MPI_STATUS_IGNORE);

    MPI_Cart_shift(comm, 0, -coord[1], &source, &dest);
    MPI_Sendrecv_replace(sb, max_a_col_b_row * local_cols, mpitype, dest, 0, source, 0, comm, MPI_STATUS_IGNORE);

    /*if (id == 3) {*/
        /*printf("coord:[%d][%d] source: %d dest: %d\n", coord[0], coord[1], source, dest);*/
        /*int local_rows = BLOCK_SIZE(coord[0], dim[0], ma);*/
        /*int local_cols = BLOCK_SIZE(coord[1], dim[1], na);*/

        /*printf("A[i][j]: \n");*/
        /*for (int i = 0; i < local_rows; ++i) {*/
            /*for (int j = 0; j < local_cols; ++j)*/
                /*printf("%.2lf ", a[i][j]);*/
            /*printf("\n");*/
        /*}*/

        /*printf("B[i][j]: \n");*/
        /*for (int i = 0; i < local_rows; ++i) {*/
            /*for (int j = 0; j < local_cols; ++j)*/
                /*printf("%.2lf ", b[i][j]);*/
            /*printf("\n");*/
        /*}*/
    /*}*/

    for (int i = 0; i < p_sqrt; ++i) {
        c = init_local_c(local_rows, local_cols, c, sc);

        /*if (id == 0) {*/

            /*printf("A[%d][%d]: \n", coord[0], coord[1]);*/
            /*for (int i = 0; i < local_rows; ++i) {*/
                /*for (int j = 0; j < local_a_col_b_row; ++j)*/
                    /*printf("%.2lf ", a[i][j]);*/
                /*printf("\n");*/
            /*}*/

            /*printf("B[%d][%d]: \n", coord[0], coord[1]);*/
            /*for (int i = 0; i < local_a_col_b_row; ++i) {*/
                /*for (int j = 0; j < local_cols; ++j)*/
                    /*printf("%.2lf ", b[i][j]);*/
                /*printf("\n");*/
            /*}*/
        /*}*/

        matmul(a, b, c, local_rows, local_a_col_b_row, local_cols);

        sum(local_rows, local_cols, c, partial_c_matrix);

        MPI_Cart_shift(comm, 1, 1, &source, &dest);
        MPI_Sendrecv_replace(sa, local_rows * max_a_col_b_row, mpitype, dest, 0, source, 0, comm, MPI_STATUS_IGNORE);

        MPI_Sendrecv_replace(&local_a_col_b_row, 1, MPI_INT, dest, 0, source, 0, comm, MPI_STATUS_IGNORE);

        MPI_Cart_shift(comm, 0, 1, &source, &dest);
        MPI_Sendrecv_replace(sb, max_a_col_b_row * local_cols, mpitype, dest, 0, source, 0, comm, MPI_STATUS_IGNORE);

    }

    /*if (id == 0) {*/

        /*printf("C[%d][%d]: \n", coord[0], coord[1]);*/
        /*for (int i = 0; i < local_rows; ++i) {*/
            /*for (int j = 0; j < local_cols; ++j)*/
                /*printf("%.2lf ", partial_c_matrix[i][j]);*/
            /*printf("\n");*/
        /*}*/
    /*}*/

    /* write the submatrix of C managed by this process */
    write_checkerboard_matrix(argv[3], (void **) partial_c_matrix, mpitype, ma, nb, comm);

    /* final timing */
    elapsed_time += MPI_Wtime();

    if (!id) {
        printf("elapsed time: %lfs\n", elapsed_time);
    }

    MPI_Finalize();

    return 0;
}






