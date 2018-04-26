// INCLUDES --------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include "papi.h"
 
// DEFINES ---------------------------------------------------------------------
#define MAXDIM 10000
#define MAXEVENTS 20
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d\n%s: line %d\n", retval,__FILE__,__LINE__); exit(retval); }
 
// VARIABLES -------------------------------------------------------------------
int a[MAXDIM][MAXDIM], b[MAXDIM][MAXDIM];
int c[MAXDIM][MAXDIM];
int dim;
int number_of_threads;

using namespace std;

// MATRIX MULTIPLICATION FUNCTIONS ---------------------------------------------
void multiplyNaive() {
        for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                        for (int k = 0; k < dim; k++) {
                                c[i][j] += a[i][k] * b[k][j];
			}
                }
        }
}

void multiplyOptimized() {
        for (int i = 0; i < dim; i++) {
		for (int k = 0; k < dim; k++) {
                	for (int j = 0; j < dim; j++) {
                                c[i][j] += a[i][k] * b[k][j];
			}
                }
        }
}

void multiplyNaiveParallel() {
	int i,j,k;

#pragma omp parallel for private(j, k) num_threads(number_of_threads)
        for (i = 0; i < dim; i++) {
                for (j = 0; j < dim; j++) {
                        for (k = 0; k < dim; k++) {
                                c[i][j] += a[i][k] * b[k][j];
			}
                }
        }
}

void multiplyOptimizedParallel() {
	int i,j,k;

#pragma omp parallel for private(j, k) num_threads(number_of_threads)
        for (i = 0; i < dim; i++) {
		for (k = 0; k < dim; k++) {
                	for (j = 0; j < dim; j++) {
                                c[i][j] += a[i][k] * b[k][j];
			}
                }
        }
}

// MATRIX FUNCTIONS ------------------------------------------------------------
void print_matrix(int matrix[][MAXDIM]) {

        for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++)
                        printf("%7d ", matrix[i][j]);
                printf("\n");
        }
        printf("\n");
}
 
void random_matrix(int matrix[][MAXDIM]) {
        for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                        matrix[i][j] = rand() % 100;
		}
	}
}

void initialize_matrix(int matrix[][MAXDIM]) {
        for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                        matrix[i][j] = 0;
		}
	}
}

// PAPI FUNCTIONS ---------------------------------------------------------------
int setupPAPI(int EventSet) {
  	int retval;

	if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT )
      		ERROR_RETURN(retval);
	
	if ((retval = PAPI_create_eventset(&EventSet)) != PAPI_OK)
      		ERROR_RETURN(retval);

   	if ((retval = PAPI_add_event(EventSet, PAPI_TOT_INS)) != PAPI_OK)
      		ERROR_RETURN(retval);

   	if ((retval = PAPI_add_event(EventSet, PAPI_TOT_CYC)) != PAPI_OK)
      		ERROR_RETURN(retval);

	if ((retval = PAPI_add_event(EventSet, PAPI_L1_TCM)) != PAPI_OK)
      		ERROR_RETURN(retval);

	if ((retval = PAPI_add_event(EventSet, PAPI_L2_TCM)) != PAPI_OK)
      		ERROR_RETURN(retval);

	if ((retval = PAPI_start(EventSet)) != PAPI_OK)
		ERROR_RETURN(retval);
}

void closePAPI(int EventSet, double real_time) {
	int retval;
	double ipc;
  	long long ins, total_flop, flops;
	long long values[MAXEVENTS];

	/* read the counter values and store them in the values array */
   	if ((retval=PAPI_read(EventSet, values)) != PAPI_OK)
      		ERROR_RETURN(retval);

   	/* Stop counting and store the values into the array */
   	if ((retval = PAPI_stop(EventSet, values)) != PAPI_OK)
      		ERROR_RETURN(retval);
	
	ins = values[0];
	ipc = ((double)values[0])/((double)values[1]);
	total_flop = 3*pow(dim,3) + dim*(dim+1);
	flops = total_flop/real_time;

	cout << "Real time:\t\t" << real_time << " seconds" << endl;
	cout << "Total instructions:\t" << ins << endl;
	cout << "Instructions per cycle:\t" << ipc << endl;
	cout << "Total FLOP:\t\t" << total_flop << endl;
	cout << "FLOP/s:\t\t\t" << flops << endl;
	cout << "L1 cache misses:\t" << values[2] << endl;
	cout << "L2 cache misses:\t" << values[3] << endl;
	cout << __FILE__ << "\t\tPASSED\n" << endl;
  	PAPI_shutdown();
}

float checkPerformance(bool regular)
{
	int EventSet = PAPI_NULL;	
	double real_time;	
	struct timespec now, tmstart;
	random_matrix(a);
        random_matrix(b);
	initialize_matrix(c);

	EventSet = setupPAPI(EventSet);

    	clock_gettime(CLOCK_MONOTONIC, &tmstart);

	if (regular) {
		if (number_of_threads > 1) {
			multiplyNaiveParallel();
		} else multiplyNaive();
	} else {
		if (number_of_threads > 1) {
			multiplyOptimizedParallel();
		} else multiplyOptimized();
	}
	clock_gettime(CLOCK_MONOTONIC, &now);
	real_time = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));

	closePAPI(EventSet, real_time);

	return real_time;
}

void multiplyMatricesNxN(int operation) {
	cout << "MULTIPLICATION OF MATRICES " << dim << "X" << dim << endl;
	float firstElapsedTime, secondElapsedTime, ratio;

	switch(operation) {
		case 1:
			cout << "\tOp1 : NxN matrix multiplication (normal) - " << number_of_threads << " threads" << endl;
			checkPerformance(true);
			break;
		case 2:
			cout << "\tOp2 : NxN matrix multiplication (alternative) - " << number_of_threads << " threads" << endl;
			checkPerformance(false);
			break;
		case 3:
			cout << "\tOp3 : ratio of performance between different algorithms" << endl;
			firstElapsedTime = checkPerformance(true);
			secondElapsedTime = checkPerformance(false);

			ratio = (float)(firstElapsedTime)/(float)(secondElapsedTime);
			cout << "Ratio between first and second cycle's execution times: " << ratio << endl << endl;
			break;		
		default:
			break;
	}

	return;
}
 
// MAIN ------------------------------------------------------------------------
int main(int argc, char* argv[]) {
	srand(time(NULL));
	size_t operation;	
	if (argc == 4) {
		istringstream ss1(argv[1]);
		istringstream ss2(argv[2]);
		istringstream ss3(argv[3]);
		
		if (!(ss1 >> operation) || !(ss2 >> number_of_threads) || !(ss3 >> dim)) {
			cerr << "Error! Arguments must be integers." << endl;
		}

		multiplyMatricesNxN(operation);
	} else {
		cout << "usage: " << argv[0] << " <operation> <number_of_threads> <size>" << endl;
		cout << "\tOp1 : NxN matrix multiplication (normal)" << endl;
		cout << "\tOp2 : NxN matrix multiplication (alternative)" << endl;
		cout << "\tOp2 : ratio of performance between different algorithms" << endl;
		cout << "\tnumber_of_threads : number of threads" << endl;
		cout << "\tsize : matrix size" << endl;
	}

	return 0;
}
