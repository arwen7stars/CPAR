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

    unsigned long size_per_process = number_odd_elements/size;
    unsigned long blockLow = (rank*number_odd_elements/size)*2 + 1;
    unsigned long blockHigh = (((rank+1)*number_odd_elements/size) - 1) * 2 + 1;
	  unsigned long sqrt_last_number = sqrt(last_number);

    std::vector<bool> sieved_vector(size_per_process, false);   // creates a vector for each process
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
      start = MPI_Wtime();
    }

    size_t startBlock;
    for (size_t k = 3; k <= static_cast<size_t>(sqrt(last_number));) {
      if (k * k < blockLow) {
        startBlock = blockLow;
        if (blockLow % k != 0) {
          startBlock += -(blockLow % k) + k;
          if (startBlock % 2 == 0) {
            startBlock += k;
          }
        }
      } else {
        startBlock = k * k;
      }

      for (size_t multiple = startBlock; multiple <= blockHigh; multiple += 2 * k) {
        sieved_vector[(multiple - blockLow) / 2] = true;
      }

      if (rank == 0) {
        k += 2;
      }
      MPI_Bcast(&k, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    }

    size_t blockPrimes = 0;
    size_t start_value = blockLow;
	  std::ofstream fout;
	  fout.open("primes.csv");
    if (rank == 0) {
      fout << "2, " << std::endl;
      start_value += 2;
      blockPrimes += 1;
    }
    for (size_t number = start_value; number <= blockHigh; number += 2) {
      if (!sieved_vector[(number - blockLow) / 2]) {
        fout << number << ", ";
        blockPrimes++;
      }
    }

    size_t AllBlocksPrimes = 0;

    MPI_Reduce(&blockPrimes, &AllBlocksPrimes, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      end = MPI_Wtime();
      std::cout << "Primes found: " << AllBlocksPrimes << std::endl;
      std::cout << (end - start) << std::endl;
    }

    fout.close();
    MPI_Finalize();

  } else {
      std::cout << "usage: "<< argv[0] << " <last_number>" << std::endl;
  }

	return 0;
}