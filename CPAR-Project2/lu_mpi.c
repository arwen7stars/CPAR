#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define MAXDIM 6000

// MATRIX FUNCTIONS ------------------------------------------------------------
void print_matrix(float *matrix, size_t no_cells, size_t dim)
{
    int i;
    for (i = 0; i < no_cells; i++) {
        printf("% 4.1f\t", matrix[i]);
        
        if ((i + 1) % dim == 0) {
            printf("\n");
        }
    }
}

float* initialize_matrix(size_t no_cells, size_t dim, int random)
{
    int i;
    float *matrix = malloc(sizeof(float) * no_cells);

    for(i = 0; i < no_cells; i++) {
        if(random) {
            matrix[i] = rand() % 100 + 1;
        } else {
            matrix[i] = 0;
        }
    }

    return matrix;
}


/*
    Auxiliar function to read_matrix_file
*/
char* remove_commas(char* str) {
    char *r, *w;
    for (w = r = str; *r; r++) {
        if (*r != ',') {
            *w++ = *r;
        }
    }
    *w = '\0';

    return str;
}

/*
    Reads the matrix that will be used in the LU decomposition algorithm
*/
float** read_matrix_file(char* filename, int* dim) {
    
    FILE *file = fopen(filename, "r");
    
    if (file)
    {
        unsigned long i, j, k;
        float* array = initialize_matrix(MAXDIM*MAXDIM, MAXDIM, 0);
        char buffer[MAXDIM*1000], *ptr;
        
        int tmp_dim = -1;
        for (i = 0; fgets(buffer, sizeof buffer, file); i++)
        {
            char* tmp = remove_commas(buffer);

            char *p = buffer;
            int j = 0;
            while (*p) {
                if (isdigit(*p)) {                              // Upon finding a digit, ...
                    if(i == 0) {
                        array[j] = (float)strtol(p, &p, 10);     // Read a number, ...
                    } else {
                        array[i*tmp_dim + j] = (float)strtol(p, &p, 10);
                    }
                    
                    if(j > tmp_dim) {
                        tmp_dim = j;
                    }

                    j++;
                } else {                                        // Otherwise, move on to the next character.
                    p++;
                }
            }
        }
        fclose(file);

        array = realloc(array, sizeof(float) * i*i);
        *dim = i;

        return array;
    }
    else
    {
        perror(filename);
	return NULL;
    }
}

/*
    Writes matrices into csv files. If integer is set to true, this will print the integer version of the matrix
*/
void write_matrix(float* matrix, int no_cells, int dim, char* filename, int integer) {
    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < no_cells; i++) {
        if(integer) {
            int tmp = (int)matrix[i];
            fprintf(f, "%d, ", tmp);
        } else {
            fprintf(f, "%f, ", matrix[i]);
        }

        if ((i+1) % dim == 0) {
            fprintf(f, "\n");
        }
    }
}

/*
    Obtains LU so that we can check whether A=LU
*/
void checkResult(size_t no_cells, size_t dim, float* L, float* U) {
    float *result = initialize_matrix(no_cells, dim, 0);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                    result[i*dim + j] += L[i*dim + k] * U[k*dim + j];
            }
        }
    }

    print_matrix(result, no_cells, dim);
    free(result);
}

/*
    MPI_BCAST implementation
*/
void my_bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator) {
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);

    if (world_rank == root) {
        // If we are the root process, send our data to everyone
        int i;
        for (i = 0; i < world_size; i++) {
            if (i != world_rank) {
                MPI_Send(data, count, datatype, i, 0, communicator);
                //printf("Sent %f to process %d from process %d.\n", data[0], i,world_rank);
            }
        }
    } else {
        // If we are a receiver process, receive the data from the root
        MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
        //printf("Received %f from process %d at process %d.\n",data[0], root,world_rank);
    }
}

/*
    Given the matrix obtained by LU decomposition, this function sets the L and U arrays properly
*/
void setLU(float* M, float* L, float* U, int dim)
{
    int i, j;

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (j > i) {
                L[i*dim + j] = 0;
                U[i*dim + j] = M[i*dim + j];
            } else if (i == j) {
                L[i*dim + j] = 1;
                U[i*dim + j] = M[i*dim + j];
            } else {
                L[i*dim + j] = M[i * dim + j];
                U[i*dim + j] = 0;
            }
        }
    }
}

void calculateRowLU(float **row, float *pivot_row, size_t dim)
{
    float factor = **row / pivot_row[0];

    int i;
    for (i = 1; i < dim; i++) {                             // goes through all coefficients in current row that are above the main diagonal...
        (*row)[i] = (*row)[i] - factor * pivot_row[i];      // ...and subtracts them by the corresponding coefficient in the pivot row multiplied by factor
    }
    **row = factor;                                         // instead of setting this coefficient to zero like we would do if were using the usual gaussian method we set this coefficient to factor
    // NOTE: **row points to A[i][main_elem], being i the same as in luDecomposition
}

int main(int argc, char *argv[])
{
   if (argc == 2) {
       char *matrix_dim = argv[1];
       char *ptr;
       
       srand(time(NULL));
       
       int dim = strtol(matrix_dim, &ptr, 10);

       if (*ptr){
            printf("Conversion error, non-convertible part: %s", ptr);
        }

        float *A = initialize_matrix(dim * dim, dim, 1);
        float *L = initialize_matrix(dim * dim, dim, 0);
        float *U = initialize_matrix(dim * dim, dim, 0);

        const int root_process = 0;
        int size, rank;
        
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        if (rank == root_process) {
            /*printf("[A]\n");
            print_matrix(A, dim * dim, dim);*/
        }

        int i;
        int main_elem;
        double start = MPI_Wtime();
        
        for (main_elem = 0; main_elem < dim-1; main_elem++) {
            float *pivot_row = &A[main_elem * dim + main_elem];             // points to element of pivot row that belongs to the main diagonal
            
            for (i = main_elem + 1; i < dim; i++) {                         // goes through all coefficients of pivot column below A[i][i]
                if (i % size == rank) {                                     // guarantees that the work is properly divided by all processes
                    float *row = &A[i * dim + main_elem];                   // points to coefficient below A[i][i]...
                    // this coefficient will be used in calculating the factor by which to subtract the current row coefficient so that it turns to zero
                    
                    calculateRowLU(&row, pivot_row, dim - main_elem);
                    //printf("data: %f\n", row[0]);
                }
            }

            for (i = main_elem + 1; i < dim; i++) {
                float *row = &A[i * dim + main_elem];                                   // pointer to row to be sent to all processes or to be received
                my_bcast(row, dim - main_elem, MPI_FLOAT, i % size, MPI_COMM_WORLD);    // keeps data synchronized in all processes
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        double end = MPI_Wtime();
        if(rank == root_process) {
            setLU(A, L, U, dim);
            /*printf("\n[L]\n");
            print_matrix(L, dim * dim, dim);

            printf("\n[U]\n");
            print_matrix(U, dim * dim, dim);

            printf("\n[LU]\n");
            checkResult(dim*dim, dim, L, U);*/

            printf("mpi: %f s\n", end - start);
        }

        free(A);
        free(L);
        free(U);
        MPI_Finalize();
        return EXIT_SUCCESS;
   } else {
       printf("usage: %s <size>\n", argv[0]);
       return EXIT_FAILURE;
   }
}
