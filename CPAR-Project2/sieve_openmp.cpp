#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <omp.h>

void sieve(unsigned long long last_number) {
    bool* prime_numbers = new bool[last_number/2];      // all elements are set to false
    size_t sqrt_last_number = sqrt(last_number);

    #pragma omp parallel for num_threads(2) schedule(dynamic)
	for (size_t i = 3; i <= sqrt_last_number; i+=2)
	{
		if (!prime_numbers[i/2])
		{
			//std::cout << std::endl << i << " " << i/2 << std::endl << std::endl;
            for (size_t j = i*i; j < last_number; j += 2*i)
			{
				//std::cout << "current multiple: " << j << " " << j/2 << " " << j/2/2 << std::endl;
				prime_numbers[j/2] = true;
			}
		}
	}

    size_t counter = 1;
    //cout << 2 << ", ";
	for(size_t i = 3; i < last_number; i+=2) {
		if(!prime_numbers[i/2]) {
			counter++;
			//std::cout << i << ", ";
		}
	}

	std::cout << "Prime numbers found: " << counter << std::endl;
    delete[] prime_numbers;
}

int main(int argc, char* argv[]) {
    if (argc == 2) {
        char *matrix_dim = argv[1];
        char *ptr;
		struct timespec now, tmstart;

        srand(time(NULL));

        long long last_number = strtoull(matrix_dim, &ptr, 10);
		
		clock_gettime(CLOCK_MONOTONIC, &tmstart);
		sieve(last_number);
		clock_gettime(CLOCK_MONOTONIC, &now);
		float real_time = (double)((now.tv_sec+now.tv_nsec*1e-9) - (double)(tmstart.tv_sec+tmstart.tv_nsec*1e-9));

		std::cout << "Elapsed time: "<< real_time << "s" << std::endl;
    } else {
        std::cout << "usage: "<< argv[0] << " <size>" << std::endl;
    }
}