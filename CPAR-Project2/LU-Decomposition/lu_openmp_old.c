#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define MAXDIM 10000

float a[MAXDIM][MAXDIM], l[MAXDIM][MAXDIM], u[MAXDIM][MAXDIM], result[MAXDIM][MAXDIM];
int dim;

// MATRIX FUNCTIONS ------------------------------------------------------------
void print_matrix(float matrix[][MAXDIM]) {
    for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                    printf("%5.0f ", matrix[i][j]);
            printf("\n");
    }
    printf("\n");
}
 
void random_matrix(float matrix[][MAXDIM]) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            matrix[i][j] = rand() % 100;
		}
	}
}

void initialize_matrix(float matrix[][MAXDIM]) {
    for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                    matrix[i][j] = 0;
    }
}
}

void checkResult() {
    for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                    for (int k = 0; k < dim; k++) {
                            result[i][j] += l[i][k] * u[k][j];
        }
            }
    }
    print_matrix(result);

}

// Uses Doolittle algorithm to decompose A matrix into lower (L) and upper (U) matrices
void lu_decomposition()
{
    
    int i,j,k;
   // #pragma omp parallel for private(j, k) num_threads(4) schedule(dynamic)
    for(i = 0; i < dim; i++) {
        l[i][i] = 1;

        //#pragma omp parallel private(k) num_threads(4)
        #pragma omp parallel for private(k) num_threads(4) schedule(dynamic)
        for(j = i; j < dim;j++) {
            float sum = 0.0;
            //#pragma omp simd
            for(k = 0;k <= i-1; k++) {
                sum += l[i][k]*u[k][j];
            }
            u[i][j] = a[i][j] - sum;
        }
    
        //#pragma omp parallel private(k) num_threads(4)
        #pragma omp parallel for private(k) num_threads(4) schedule(dynamic)
        for(j=i+1; j<dim; j++) {
            float sum = 0.0;
            //#pragma omp parallel for
              //#pragma omp simd

            for(k=0; k <= i-1;k++) {
                sum += l[j][k]*u[k][i];
            }
            l[j][i] = (a[j][i]-sum) / u[i][i];
        }
    }
}

long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC;
    return elapsed;
}

int main(int argc, char* argv[]) {
    if (argc == 2) {
        char *matrix_dim = argv[1];
        char *ptr;
        struct timespec now, tmstart;

        srand(time(NULL));

        dim = strtol(matrix_dim, &ptr, 10);

        random_matrix(a);
        initialize_matrix(l);
        initialize_matrix(u);

        clock_gettime(CLOCK_MONOTONIC, &tmstart);
        lu_decomposition();
        clock_gettime(CLOCK_MONOTONIC, &now);

        float real_time = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));

        printf("Elapsed time: %f s\n", real_time);

        //print_matrix(a);

        //initialize_matrix(result);
        //checkResult();
    } else {
        printf("usage: %s <size>\n", argv[0]);
    }
}
