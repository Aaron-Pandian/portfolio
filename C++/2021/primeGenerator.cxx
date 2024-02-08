// Aaron Pandian (anp3238)
#include <iostream>
#include <bits/stdc++.h>
#include "primegenerator.cxx"
#include <vector>
using std::vector;
using namespace std;

int main()
{
    // call number of even numbers to test
    int even_numbers;
    cin >> even_numbers;
    // loop over the even numbers
    for (int i = 4; i <= (even_numbers * 2) + 2; i += 2)
    {
        primegenerator sequence;
        vector<int> primes;
        primes.push_back(sequence.nextprime());
        int counter = 0;
        // for even number, generate all primes, I am storing all values in array
        while (primes[counter] <= i)
        {
            int number = sequence.nextprime();
            primes.push_back(number);
            counter++;
        }
        // create a test to test goldbach conjecture against prime generator
        for (int j = 0; primes[j] <= i; j++)
        {
            int diff = i - primes[j];
            // Check if the difference is also a prime number
            if (is_prime(diff))
            {
                // print triple of q p r
                int r = (primes[j] + diff) / 2;
                cout << primes[j] << "," << diff << "," << r << endl;
                /* old print
        cout << "The number " << i << " is " << primes[j] << "+" << diff << endl;
         */
                break;
            }
        }
    }
    return 0;
}
