#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

#include <mpi.h>

int main(int argc, char* argv[]) {
  if (argc == 2) {
    char *matrix_dim = argv[1];
    char *ptr;
		double start = 0, end = 0;

    srand(time(NULL));

    long long last_number = strtoull(matrix_dim, &ptr, 10);
    int rank, size;

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    size_t number_odd_elements = last_number / 2;

	  unsigned long sqrt_last_number = sqrt(last_number);
    unsigned long size_per_process = number_odd_elements/size;              // divides the array by all processes

    unsigned long lowerBound = (rank*size_per_process)*2 + 1;               // lower bound of the array of proccess with id "rank"
    unsigned long upperBound = ((rank+1)*number_odd_elements)/size*2 - 1;   // upper bound of the array of proccess with id "rank"
    // note that elements with index i/2 have an value of [i] being that the reason why both lower and upper bounds are multiplied by 2

    std::vector<bool> sieved_vector(size_per_process, false);               // creates a vector for each process
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
      start = MPI_Wtime();
    }

    size_t startBlock;
    for (size_t prime = 3; prime <= sqrt_last_number; ) {
      if (prime * prime < lowerBound) {
        startBlock = lowerBound;
        if (lowerBound % prime != 0) {
          startBlock += -(lowerBound % prime) + prime;
          if (startBlock % 2 == 0) {
            startBlock += prime;
          }
        }
      } else {
        startBlock = prime * prime;
      }

      for (size_t multiple = startBlock; multiple <= upperBound; multiple += 2 * prime) {
        sieved_vector[(multiple - lowerBound) / 2] = true;
      }

      if (rank == 0) {
        prime += 2;
      }

      MPI_Bcast(&prime, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    size_t blockPrimes = 0;
    size_t start_value = lowerBound;

    if (rank == 0) {
      std::cout << "2, " << std::endl;
      start_value += 2;
      blockPrimes += 1;
    }

    for (size_t number = start_value; number <= upperBound; number += 2) {
      if (!sieved_vector[(number - lowerBound) / 2]) {
        std::cout << number << ", ";
        blockPrimes++;
      }
    }

    size_t totalPrimes = 0;

    MPI_Reduce(&blockPrimes, &totalPrimes, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      end = MPI_Wtime();
      
      std::cout << "Primes found: " << totalPrimes << std::endl;
      std::cout << "Elapsed time: " << end - start << "s" << std::endl;
    }

    MPI_Finalize();

  } else {
      std::cout << "usage: "<< argv[0] << " <last_number>" << std::endl;
  }

	return 0;
}