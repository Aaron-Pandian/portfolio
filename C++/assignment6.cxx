#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>
using std::abs;
using std::pow;
using std::sqrt;
#include <functional>
using std::function;

// Create lambda funciton to apprixmiate the gradient (using equation fiven in lecture) of an input double repres\
enting a function w newton method using h=1e-6 and then pass it through the find_zero function.

double find_zero(function<double(double)> f, function<double(double)> g)
{
    double x = sqrt(f(x) + 2);
    //  while (abs(f(x)/fprime(x)) >= .1) {
    for (int i = x; i < 5; i++)
    {
        /*
        double h = pow(10.0, -6);
        double g = (f(x+h)-f(x))/h;
        */
        x = x - f(x) / g(x); // applying fprime, so get rid of 'g'
        cout << x << endl;
    }
    return x;
}

int main()
{

    double x;
    cout << "Enter a number for x: ";
    cin >> x;

    auto sqrt2 = find_zero([x](double arb) -> double
                           { return x * x - 2; },
                           [x](double arb) -> double
                           { return 2 * x; });

    cout << "The root of the number is " << sqrt2 << '\n';

    return 0;
}