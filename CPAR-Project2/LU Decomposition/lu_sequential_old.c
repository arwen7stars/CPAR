#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXDIM 10000

float a[MAXDIM][MAXDIM], l[MAXDIM][MAXDIM], u[MAXDIM][MAXDIM], result[MAXDIM][MAXDIM];
int dim;

// MATRIX FUNCTIONS ------------------------------------------------------------
void print_matrix(float matrix[][MAXDIM]) {
    for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                    printf("%5.1f ", matrix[i][j]);
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
    for(int i = 0; i < dim; i++) {
        l[i][i] = 1;
 
        for(int j = i; j < dim;j++) {
            float sum = 0.0;

            for(int k = 0;k <= i-1; k++) {
                sum += l[i][k]*u[k][j];
            }
            u[i][j] = a[i][j] - sum;
        }
    
        for(int j=i+1; j<dim; j++) {
            float sum = 0.0;
            
            for(int k=0; k <= i-1;k++) {
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
        clock_t t1, t2, elapsed;

        srand(time(NULL));

        dim = strtol(matrix_dim, &ptr, 10);

        random_matrix(a);
        initialize_matrix(l);
        initialize_matrix(u);

        t1 = clock();
        lu_decomposition();
        t2 = clock();

        elapsed = timediff(t1, t2);

        printf("Elapsed time: %ld s\n", elapsed);


        /*print_matrix(a);
        print_matrix(l);
        print_matrix(u);

        initialize_matrix(result);
        checkResult();*/
    } else {
        printf("usage: %s <size>\n", argv[0]);
    }
}
