#include "Polynomial.h"
#include <iostream>
#include <vector>

using namespace std;

int main() {
    // Input for Polynomial 1
    int degree1;
    cout << "Enter the degree and coefficients for Polynomial 1 (first the degree, then the coefficients):\n";
    cin >> degree1;

    vector<double> coeffs1(degree1 + 1); // Vector to store coefficients
    for (int i = 0; i <= degree1; ++i) {
        cin >> coeffs1[i];
    }

    // Create Polynomial 1
    Polynomial p1(coeffs1);

    // Input for Polynomial 2
    int degree2;
    cout << "Enter the degree and coefficients for Polynomial 2 (first the degree, then the coefficients):\n";
    cin >> degree2;

    vector<double> coeffs2(degree2 + 1); // Vector to store coefficients
    for (int i = 0; i <= degree2; ++i) {
        cin >> coeffs2[i];
    }

    // Create Polynomial 2
    Polynomial p2(coeffs2);

    // Print Polynomials
    cout << "Polynomial 1: " << p1 << endl;
    cout << "Polynomial 2: " << p2 << endl;

    // Sum
    Polynomial sum = p1 + p2;
    cout << "Sum: " << sum << endl;

    // Subtraction
    Polynomial difference = p1 - p2;
    cout << "Subtraction: " << difference << endl;

    // Product
    Polynomial product = p1 * p2;
    cout << "Product: " << product << endl;

    // Degree of Polynomial 1
    cout << "Degree of Polynomial: " << p1.degree() << endl;

    // Evaluation of Polynomial 1 at x = 2
    double x = 2;
    cout << "Evaluation of polynomial at x=" << x << ": " << p1.evaluate(x) << endl;

    // Derivative of Polynomial 1
    Polynomial derivative = p1.derivative();
    cout << "Derivative: " << derivative << endl;

    // Integral of Polynomial 1
    Polynomial integral = p1.integral();
    cout << "Integral: " << integral << endl;

    // Definite integration from x = 0 to x = 1
    cout << "Definite integration from x=0 to x=1: " << p1.integral(0, 1) << endl;

    // Real roots of Polynomial 1 (approximately, with Newton's method)
    cout << "Real roots of Polynomial 1: [" << p1.getRoot(-2) << ", " << p1.getRoot(1) << "]" << endl;

    return 0;
}
