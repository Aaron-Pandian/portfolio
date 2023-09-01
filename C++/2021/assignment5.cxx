#include <iostream>
#include <vector>
using std::cin;
using std::cout;
using std::vector;
using namespace std;

/* Write class pascal with method pascal(n) will output a pascal triangle with n rows. Also include a method print(int m) that prints the triangle but replaces coefficients % m != 0 with a *, and a space if equal to zero. */

class pascal
{
    int n;
    vector<int> elements;

public:
    pascal(int rows) : n(rows), elements(vector<int>(n * n)){};

    int get(int i, int j)
    {
        int x = 1;
        // calculate values for i,j
        for (int k = 1; k <= j; k++)
        {
            x = x * (i - k + 1) / k;
        }
        return x;
    };

    // Method for printing pascals triangle
    void print()
    {
        for (int i = 1; i <= n; i++)
        {

            for (int j = 1; j < n - i + 1; j++)
                cout << "  ";

            int val = 1;
            for (int j = 1; j <= i; j++)
            {
                cout << val << "   ";
                val = val * (i - j) / j;
            }
            cout << endl;
        }
    };

    // Method for stars
    void print(int m)
    {
        if (m == 0)
        {
            // Stops Program
            return;
        }
        for (int i = 0; i < n; i++)
        {

            for (int j = 0; j < n - i - 1; j++)
            {
                cout << " ";
            }

            for (int j = 0; j <= i; j++)
            {
                int num = get(i, j);
                if (num % m == 0)
                {
                    cout << "  ";
                }
                else
                {
                    cout << " *";
                }
            }
            cout << endl;
        }
    };
};

int main()
{
    int m;
    int n;
    cout << "Enter the rows: ";
    cin >> n;
    for (int i = 0; i, 3; i++)
    {
        // Print the pascal display for given n
        pascal triangle(n);
        triangle.print();
        cout << "Enter modulo number: ";
        cin >> m;
        triangle.print(m);
        if (m == 0)
        {
            break;
        }
    }
    return 0;
}
