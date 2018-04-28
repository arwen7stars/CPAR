#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

void write_primes(std::vector<bool> prime_numbers, unsigned long long last_number) {
	std::ofstream fout;
	fout.open("primes.csv");
	size_t counter = 1;

    fout << 2 << ", ";
	for(size_t i = 3; i < last_number; i+=2) {		// even numbers are ignored (4,6,8,...)
		if(!prime_numbers[i/2]) {					// prime numbers in this array are the false ones
			counter++;
			fout << i << ", ";
		}
	}

	fout.close();

	std::cout << "Prime numbers found: " << counter << std::endl;

}

/*
	"std::vector<bool> is a possibly space-efficient specialization of std::vector for the type bool."
	http://en.cppreference.com/w/cpp/container/vector_bool
*/
void sieve(unsigned long long last_number) {
	std::vector<bool> prime_numbers(last_number/2, false);		// the amount of prime numbers is smaller than the amount of uneven numbers (since there is no prime number that is even except two)
	unsigned long sqrt_last_number = sqrt(last_number);

	// in this implementation, even numbers are ignored (note that the only number that is prime and even is 2)
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
        std::cout << "usage: "<< argv[0] << " <last_number>" << std::endl;
    }
}
