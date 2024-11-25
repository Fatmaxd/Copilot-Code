#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class Polynomial {
private:
    vector<double> coeffs; // Store coefficients of the polynomial

public:
    // Constructors
    Polynomial();
    Polynomial(const vector<double>& coefficients);
    Polynomial(const Polynomial& other);

    // Destructor
    ~Polynomial();

    // Assignment operator
    Polynomial& operator=(const Polynomial& other);

    // Arithmetic operators
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;

    // Equality operator
    bool operator==(const Polynomial& other) const;

    // Output operator
    friend ostream& operator<<(ostream& out, const Polynomial& poly);

    // Utility functions
    int degree() const;
    double evaluate(double x) const;
    Polynomial derivative() const;
    Polynomial integral() const;
    double integral(double x1, double x2) const;
    double getRoot(double guess = 1, double tolerance = 1e-6, int maxIter = 100) const;

    void setCoefficients(const vector<double>& coefficients);
    double getCoefficient(int degree) const;
};

#endif // POLYNOMIAL_H
