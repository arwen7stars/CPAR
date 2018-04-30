#define _POSIX_C_SOURCE 199309L
#define __STDC_WANT_LIB_EXT2__ 1

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <omp.h>

#define MAXDIM 6000

// MATRIX FUNCTIONS ------------------------------------------------------------
void print_matrix(float** matrix, int dim) {
    for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                    printf("%6.1f ", matrix[i][j]);
            printf("\n");
    }
    printf("\n");
}


float** initialize_matrix(int dim, int random) {
    int i, j;

    float **matrix = malloc(sizeof(float *) * dim);

    for (i = 0; i < dim; i++) {
        matrix[i] = malloc(sizeof(float) * dim);
        
        for (j = 0; j < dim; j++) {
            if (random) {               // true is 1
                matrix[i][j] = rand() % 100 + 1;
            } else {
                matrix[i][j] = 0;
            }
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
        float** array = initialize_matrix(MAXDIM, 0);
        char buffer[MAXDIM*1000], *ptr;
        
        for (i = 0; fgets(buffer, sizeof buffer, file); i++)
        {
            char* tmp = remove_commas(buffer);

            char *p = buffer;
            int j = 0;
            while (*p) {
                if (isdigit(*p)) {                              // Upon finding a digit, ...
                    array[i][j] = (float)strtol(p, &p, 10);     // Read a number, ...
                    j++;
                } else {                                        // Otherwise, move on to the next character.
                    p++;
                }
            }
        }
        fclose(file);

        array = realloc(array, sizeof(float *) * i);
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
void write_matrix(float** matrix, int dim, char* filename, int integer) {
    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        return;
    }

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++){
            if(integer) {
                int tmp = (int)matrix[i][j];
                if(j-1 == dim) {
                    fprintf(f, "%d", tmp);
                } else {
                    fprintf(f, "%d, ", tmp);
                }
            } else {
                if(j-1 == dim) {
                    fprintf(f, "%f", matrix[i][j]);
                } else {
                    fprintf(f, "%f, ", matrix[i][j]);
                }
            }
        }
        fprintf(f, "\n");
    }
}

/*
    Obtains LU so that we can check whether A=LU
*/
void checkResult(float** l, float** u, int dim) {
    float** result = initialize_matrix(dim, 0);

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
            result[i][j] += l[i][k] * u[k][j];
            }
        }
    }

    print_matrix(result, dim);
}

/*
    Given the matrix obtained by lu_decomposition, this function sets the L and U arrays properly
*/
void setLU(float** matrix, float** l, float** u, int dim, int no_threads) {
    
    int i, j;
    
    #pragma omp parallel for private(j) num_threads(4) schedule(dynamic)
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            if (i < j) {
                l[i][j] = 0;
                u[i][j] =  matrix[i][j];
            } else if (i > j) {
                l[i][j] = matrix[i][j];
                u[i][j] = 0;
            } else {
                l[i][j] = 1;
                u[i][j] = matrix[i][j];
            }
      }
   }
}

/*
    LU-decomposition based on Gaussian Elimination
    Obtains the lower and upper triangular matrices only by using the original matrix
    http://www.personal.psu.edu/jhm/f90/lectures/lu.html

    About the function lu_decomposition...
    Arguments: a -> original square matrix; dim -> matrix size
    Returns new square matrix which be very easily decomposed into L and U
*/
int lu_decomposition(float** a, int dim, int no_threads)
{
    int i, j, k;
    int success = 1;

    for (i = 0; i < dim; i++) {
        #pragma omp parallel for num_threads(no_threads) schedule(dynamic)
        for (j = i + 1; j < dim; j++) {
            float factor = a[j][i] / a[i][i];   // factor by which to multiply the equation of row i so that a[j][i] becomes zero
            
            for (k = i + 1; k < dim; k++) {     // goes through all coefficients of row j that are above the main diagonal and substracts them by the corresponding coefficient in row i
            // Note: each coefficient corresponds to column k
                a[j][k] -= factor * a[i][k];
            }
            a[j][i] = factor;                   // instead of setting a[j][i] to zero like we do usually in gaussian elimination we keep this value so that we obtain the lower triangular
        }
    }

    return success;
}

int main(int argc, char* argv[]) {
    if (argc == 3) {
        srand(time(NULL));
        struct timespec now, tmstart;

        /*float** A = initialize_matrix(matrix_size, 1);
        write_matrix(A, matrix_size, "example.csv", 1);*/

        /*char *matrix_dim = argv[1];
	    char *nthreads = argv[2];
        char *ptr_dim;
	    char *ptr_nt;

        int dim = strtol(matrix_dim, &ptr_dim, 10);
	    int no_threads = strtol(nthreads, &ptr_nt, 10);

        if (*ptr_dim) {
            printf("ERROR: Conversion error, non-convertible part: %s", ptr_dim);
            return EXIT_FAILURE;
        } else if(*ptr_nt) {
            printf("ERROR: Conversion error, non-convertible part: %s", ptr_nt);
            return EXIT_FAILURE;
	    }
        
        float** A = initialize_matrix(dim, 1);*/
        
        char* filename = argv[1];
        char* nthreads = argv[2];
        char* matrix_size = argv[3];
        char *ptr;

        int dim;
        float** A = read_matrix_file(filename, &dim);

        if(A == NULL) {
            perror("Error reading matrix file.\n");
	        return EXIT_FAILURE;	
	    }

        int no_threads = strtol(nthreads, &ptr, 10);

        if (*ptr) {
            printf("ERROR: Conversion error, non-convertible part: %s", ptr);
            return EXIT_FAILURE;
        }

        float** L = initialize_matrix(dim, 0);
        float** U = initialize_matrix(dim, 0);

        /*printf("[A]\n");
        print_matrix(A, dim);*/

        clock_gettime(CLOCK_MONOTONIC, &tmstart);

        lu_decomposition(A, dim, no_threads);
        setLU(A, L, U, dim, no_threads);

        clock_gettime(CLOCK_MONOTONIC, &now);

        float real_time = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));

        write_matrix(L, dim, "L.csv", 0);
        write_matrix(U, dim, "U.csv", 0);
        
        /*printf("\n[L]\n");
        print_matrix(L, dim);
        
        printf("\n[U]\n");
        print_matrix(U, dim);

        printf("\n[LU]\n");
        checkResult(L, U, dim);*/

        printf("Elapsed time: %f s\n", real_time);
    } else {
        printf("usage: %s <filename> <no_threads>\n", argv[0]);
        return EXIT_FAILURE;
    }
}
