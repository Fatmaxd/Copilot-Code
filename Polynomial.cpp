#include "Polynomial.h"

// Default constructor
Polynomial::Polynomial() : coeffs(1, 0.0) {}

// Constructor with coefficients
Polynomial::Polynomial(const vector<double>& coefficients) : coeffs(coefficients) {
    // Remove trailing zeros
    while (!coeffs.empty() && coeffs.back() == 0.0) {
        coeffs.pop_back();
    }
}

// Copy constructor
Polynomial::Polynomial(const Polynomial& other) : coeffs(other.coeffs) {}

// Destructor
Polynomial::~Polynomial() {}

// Assignment operator
Polynomial& Polynomial::operator=(const Polynomial& other) {
    if (this != &other) {
        coeffs = other.coeffs;
    }
    return *this;
}

// Addition operator
Polynomial Polynomial::operator+(const Polynomial& other) const {
    int maxDegree = max(degree(), other.degree());
    vector<double> result(maxDegree + 1, 0.0);

    for (int i = 0; i <= degree(); ++i) {
        result[i] += coeffs[i];
    }
    for (int i = 0; i <= other.degree(); ++i) {
        result[i] += other.coeffs[i];
    }

    return Polynomial(result);
}

// Subtraction operator
Polynomial Polynomial::operator-(const Polynomial& other) const {
    int maxDegree = max(degree(), other.degree());
    vector<double> result(maxDegree + 1, 0.0);

    for (int i = 0; i <= degree(); ++i) {
        result[i] += coeffs[i];
    }
    for (int i = 0; i <= other.degree(); ++i) {
        result[i] -= other.coeffs[i];
    }

    return Polynomial(result);
}

// Multiplication operator
Polynomial Polynomial::operator*(const Polynomial& other) const {
    int n = degree() + other.degree();
    vector<double> result(n + 1, 0.0);

    for (int i = 0; i <= degree(); ++i) {
        for (int j = 0; j <= other.degree(); ++j) {
            result[i + j] += coeffs[i] * other.coeffs[j];
        }
    }

    return Polynomial(result);
}

// Equality operator
bool Polynomial::operator==(const Polynomial& other) const {
    return coeffs == other.coeffs;
}

// Output operator
ostream& operator<<(ostream& out, const Polynomial& poly) {
    if (poly.coeffs.empty()) {
        out << "0";
        return out;
    }

    bool firstTerm = true;
    for (int i = poly.degree(); i >= 0; --i) {
        if (poly.coeffs[i] != 0.0) {
            if (!firstTerm) {
                out << (poly.coeffs[i] > 0 ? " + " : " - ");
            }
            firstTerm = false;
            if (abs(poly.coeffs[i]) != 1 || i == 0) {
                out << abs(poly.coeffs[i]);
            }
            if (i > 0) {
                out << "x";
                if (i > 1) {
                    out << "^" << i;
                }
            }
        }
    }

    return out;
}

// Get the degree of the polynomial
int Polynomial::degree() const {
    return coeffs.size() - 1;
}

// Evaluate the polynomial at a given value of x
double Polynomial::evaluate(double x) const {
    double result = 0.0;
    for (int i = degree(); i >= 0; --i) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Get the derivative of the polynomial
Polynomial Polynomial::derivative() const {
    if (coeffs.empty()) {
        return Polynomial();
    }

    vector<double> result(degree());
    for (int i = 1; i <= degree(); ++i) {
        result[i - 1] = coeffs[i] * i;
    }

    return Polynomial(result);
}

// Get the integral of the polynomial
Polynomial Polynomial::integral() const {
    if (coeffs.empty()) {
        return Polynomial();
    }

    vector<double> result(degree() + 2, 0.0);
    for (int i = 0; i <= degree(); ++i) {
        result[i + 1] = coeffs[i] / (i + 1);
    }

    return Polynomial(result);
}

// Get the definite integral of the polynomial from x1 to x2
double Polynomial::integral(double x1, double x2) const {
    Polynomial integralPoly = integral();
    return integralPoly.evaluate(x2) - integralPoly.evaluate(x1);
}

// Find a root of the polynomial using Newton's method
double Polynomial::getRoot(double guess, double tolerance, int maxIter) const {
    for (int i = 0; i < maxIter; ++i) {
        double f = evaluate(guess);
        double df = derivative().evaluate(guess);
        if (abs(df) < 1e-10) {
            break; // Avoid division by zero
        }
        double newGuess = guess - f / df;
        if (abs(newGuess - guess) < tolerance) {
            return newGuess;
        }
        guess = newGuess;
    }
    return guess; // Return the last guess if convergence is not achieved
}

// Set the coefficients of the polynomial
void Polynomial::setCoefficients(const vector<double>& coefficients) {
    coeffs = coefficients;
}

// Get the coefficient of a specific degree
double Polynomial::getCoefficient(int degree) const {
    if (degree < 0 || degree > this->degree()) {
        return 0.0;
    }
    return coeffs[degree];
}
