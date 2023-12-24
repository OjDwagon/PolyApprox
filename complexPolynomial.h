#ifndef COMPLEX_POLYNOMIAL_H
#define COMPLEX_POLYNOMIAL_H

#include "genericPolynomial.h"
#include "complex.h"
#include <vector>
#include <iostream>
#include <functional>

/**
 * A class the models complex valued polynomials of a non-negative degree and their operations (excluding division)
*/
class ComplexPolynomial: public GenericPolynomial<Complex, ComplexPolynomial>
{
    public:
        /**
         * Default constructor for a ComplexPolynomial, defaults to constant 0
        */
        ComplexPolynomial();

        /**
         * ComplexPolynomial constructor given coefficients
         * 
         * @param coefficients a vector of coefficeints, the ith element of the vector is used as the coefficient for x^i
        */
        ComplexPolynomial(const std::vector<double>& coefficients);

        /**
         * ComplexPolynomial constructor given coefficients
         * 
         * @param coefficients a vector of coefficeints, the ith element of the vector is used as the coefficient for x^i
        */
        ComplexPolynomial(const std::vector<Complex>& coefficients);

        /**
         * ComplexPolynomial constructor given coefficients
         * 
         * @param coefficients an initialization list of type Complex, the ith element of the list is used as the coefficient for x^i
        */
        ComplexPolynomial(const std::initializer_list<Complex>& coefficients);
        ComplexPolynomial operator*(const ComplexPolynomial& rhs) const;
        ComplexPolynomial operator*(Complex rhs) const;
        friend std::ostream& operator<<(std::ostream& out, const ComplexPolynomial& rhs);

        /**
         * Returns the product of 2 polynomials, implemented using the Discrete Fourier Transform in loglinear time
         * 
         * @param lhs the first polynomial to multiply
         * @param rhs the second polynomial to multiply
         * 
         * @returns the polynomial lhs * rhs
        */
        static ComplexPolynomial discreteFourierTransformMultiplication(const ComplexPolynomial& lhs, const ComplexPolynomial& rhs);
    protected:
        static std::string getPrintableComplex(const Complex& num);

        /**
         * Evaluates the given polynomial on the nth roots of unity, recursively using FFT
         * 
         * @param f the function to evaluate
         * @param n which roots of unity to evaluate on, located in the roots of unity vector, must be a power of 2
         * @param rootsOfUnity a shared array of roots of unity to prevent them from being recalculated
         * @returns a vector of the Complex numbers storing f evaluated at the roots of unity (in order going counter-clockwise around the unit circle)
        */
        static std::vector<Complex> evaluateAtRootsOfUnity(const ComplexPolynomial& f, unsigned int n, const std::vector<std::vector<Complex>>& rootsOfUnity); 

        /**
         * Calculates the nth roots of unity
         * 
         * @param n the number of roots of unity to calculate
         * @returns a vector of the roots of unity (ordered going counter-clockwise around the unit circle starting from 0)
        */
        static std::vector<Complex> calcRootsOfUnity(unsigned int n);

        /**
         * If n is not a power of 2, then returns n rounded up to the next power of 2
         * 
         * @param n the number to round up
         * @returns the rounded up number
        */
        static unsigned int getNextPowerOfTwo(unsigned int n);

        friend GenericPolynomial<Complex, ComplexPolynomial>;

        static std::function<bool(Complex, Complex)> comp;
};

/**
 * Adds a double to a ComplexPolynomial, the double will be treated as a degree 0 ComplexPolynomial
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator+(Complex lhs, const ComplexPolynomial& rhs);

/**
 * Subtracts a double to a ComplexPolynomial, the double will be treated as a degree 0 ComplexPolynomial
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator-(Complex lhs, const ComplexPolynomial& rhs);

/**
 * Multiplies a ComplexPolynomial by a constant
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator*(Complex lhs, const ComplexPolynomial& rhs);

/**
 * Equality check of a ComplexPolynomial with a constant, returns true if the ComplexPolynomial is degree 0 and has a matching coefficient
 * 
 * @relatesalso ComplexPolynomial
*/
bool operator==(Complex lhs, const ComplexPolynomial& rhs);

/**
 * Inequality check of a ComplexPolynomial with a constant, returns false if the ComplexPolynomial is degree 0 and has a matching coefficient
 * 
 * @relatesalso ComplexPolynomial
*/
bool operator!=(Complex lhs, const ComplexPolynomial& rhs);


#endif
