#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <omp.h>

void write_primes(bool prime_numbers[], unsigned long long last_number) {
	std::ofstream fout;
	fout.open("primes.csv");
	size_t counter = 1;

    fout << 2 << ", ";
	for(size_t i = 3; i < last_number; i+=2) {
		if(!prime_numbers[i/2]) {
			counter++;
			fout << i << ", ";
		}
	}

	fout.close();

	std::cout << "Prime numbers found: " << counter << std::endl;

}

void sieve(unsigned long long last_number, int no_threads) {
    bool* prime_numbers = new bool[last_number/2];      		// all elements are set to false
    size_t sqrt_last_number = sqrt(last_number);

    #pragma omp parallel for num_threads(no_threads) schedule(dynamic)
	for (size_t i = 3; i <= sqrt_last_number; i+=2)
	{
		if (!prime_numbers[i/2])								// every element with index i/2 corresponds in fact to the number i in the array prime_numbers
		{
            for (size_t j = i*i; j < last_number; j += 2*i)		// goes through all multiples of i that are not even numbers (e.g por i=3 => j=9,15,21...) and crosses them from the list (setting the corresponding index to true in prime_numbers)
			{
				prime_numbers[j/2] = true;
			}
		}
	}

	write_primes(prime_numbers, last_number);

    delete[] prime_numbers;
}

int main(int argc, char* argv[]) {
    if (argc == 3) {
		srand(time(NULL));

		struct timespec now, tmstart;
		char *matrix_dim = argv[1];
	    char *nthreads = argv[2];
        char *ptr_dim;
	    char *ptr_nt;

        long long last_number = strtoull(matrix_dim, &ptr_dim, 10);
	    int no_threads = strtol(nthreads, &ptr_nt, 10);

        if (*ptr_dim) {
            printf("ERROR: Conversion error, non-convertible part: %s", ptr_dim);
            return EXIT_FAILURE;
        } else if(*ptr_nt) {
            printf("ERROR: Conversion error, non-convertible part: %s", ptr_nt);
            return EXIT_FAILURE;
		}
		
		clock_gettime(CLOCK_MONOTONIC, &tmstart);
		sieve(last_number, no_threads);
		clock_gettime(CLOCK_MONOTONIC, &now);
		float real_time = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));

		std::cout << "Elapsed time: "<< real_time << "s" << std::endl;
    } else {
        std::cout << "usage: "<< argv[0] << " <last_number> <no_threads>" << std::endl;
    }
}