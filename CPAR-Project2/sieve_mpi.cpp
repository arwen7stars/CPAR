#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

#include <mpi.h>

int main(int argc, char* argv[]) {
  if (argc == 2) {
    remove("primes.csv");                 // removes "primes.csv" file if it already exists
    char *matrix_dim = argv[1];
    char *ptr;
		double start = 0, end = 0;

    srand(time(NULL));

    long long last_number = strtoull(matrix_dim, &ptr, 10);
    int rank, size;

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    size_t number_odd_elements = last_number / 2;                           // the number of primes numbers is never bigger than the number of odd numbers

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

    size_t start_i;
    for (size_t prime = 3; prime <= sqrt_last_number; ) {
      if (prime * prime < lowerBound) {       // if the multiple of the current "prime" is smaller than the lowerBound...
        start_i = lowerBound;                 // then the start value is equal to the lower bound of the array
        
        if (lowerBound % prime != 0) {              // if the lower bound of the array is not divisible by the current "prime" number
          start_i += -(lowerBound % prime) + prime; // then we make it divisible by the prime number
          
          if (start_i % 2 == 0) {             // if the start value is divisible by 2
            start_i += prime;                 // then we add the current prime number so that it isn't divisible by 2 anymore
            // note that we ignore all even numbers since we know all of them are not prime for the exception of 2 (for example if start_i was 18 before and the current prime was 3, then now it's 21)
          }
        }
      } else {
        start_i = prime * prime;  // if the multiple of the current "prime" is bigger than or equal to the lowerBound...
        // ...then the start value will be equal to the multiple of the current "prime" (sieve algorithm)
      }

      for (size_t multiple = start_i; multiple <= upperBound; multiple += 2 * prime) {    // for prime = 3, we get: multiple = 9, 15, 21, 27, etc.
        sieved_vector[(multiple - lowerBound) / 2] = true;
      }

      if (rank == 0) {
        prime += 2;
      }

      MPI_Bcast(&prime, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    size_t blockPrimes = 0;
    size_t start_value = lowerBound;

    std::fstream fout;
    fout.open ("primes.csv", std::fstream::out | std::fstream::app);

    if (rank == 0) {
      fout << "2, ";
      start_value += 2;
      blockPrimes += 1;
    }

    size_t number = start_value;
    for (size_t number = start_value; number <= upperBound; number += 2) {
      if (!sieved_vector[(number - lowerBound) / 2]) {
        fout << number << ", ";
        blockPrimes++;
      }
    }

    size_t totalPrimes = 0;

    MPI_Reduce(&blockPrimes, &totalPrimes, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      end = MPI_Wtime();
      
      std::cout << "Prime numbers found: " << totalPrimes << std::endl;
      std::cout << "Elapsed time: " << end - start << "s" << std::endl;
    }

    MPI_Finalize();

  } else {
      std::cout << "usage: "<< argv[0] << " <last_number>" << std::endl;
  }

	return 0;
}