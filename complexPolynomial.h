#ifndef COMPLEX_POLYNOMIAL_H
#define COMPLEX_POLYNOMIAL_H

#include <vector>
#include "complex.h"
#include <string>
#include <iostream>

class ComplexPolynomial 
{
    public:
        ComplexPolynomial(); // Defaults to the constant 0
        ComplexPolynomial(const std::vector<Complex>& coefficients);
        void setCoefficient(unsigned int power, Complex val); // Does nothing if an invalid power is given

        /**
         * A linear time evaluation method of the polynomial
         * 
         * @tparam V the type of the input into the polynomial, must have defined arithmetic operations with Complex Numbers
         * @param x the value that the polynomial will be evaluated at
         * @return the value of the polynomial evaluated at x as a complex number
        */
        template <typename V>
        Complex evaluate(V x) const;

        /**
         * Operator() that allows for evaluation of a polynomial in function style syntax
        */
        template <typename V>
        Complex operator()(V x) const;

        ComplexPolynomial operator+(const ComplexPolynomial& rhs) const; // Adds 2 polynomials
        ComplexPolynomial operator+(int rhs) const; // Adds a polynomial and a constant
        ComplexPolynomial operator+(double rhs) const;
        ComplexPolynomial operator+(const Complex& rhs) const;
        ComplexPolynomial& operator+=(const ComplexPolynomial& rhs);
        ComplexPolynomial& operator+=(int rhs);
        ComplexPolynomial& operator+=(double rhs);
        ComplexPolynomial& operator+=(const Complex& rhs);        
        ComplexPolynomial operator-(const ComplexPolynomial& rhs) const;
        ComplexPolynomial operator-(int rhs) const;
        ComplexPolynomial operator-(double rhs) const;
        ComplexPolynomial operator-(const Complex& rhs) const;        
        ComplexPolynomial& operator-=(const ComplexPolynomial& rhs);
        ComplexPolynomial& operator-=(int rhs);
        ComplexPolynomial& operator-=(double rhs);
        ComplexPolynomial& operator-=(const Complex& rhs);        
        ComplexPolynomial operator*(const ComplexPolynomial& rhs) const;
        ComplexPolynomial operator*(int rhs) const;
        ComplexPolynomial operator*(double rhs) const;
        ComplexPolynomial operator*(const Complex& rhs) const;        
        ComplexPolynomial& operator*=(const ComplexPolynomial& rhs);
        ComplexPolynomial& operator*=(int rhs);
        ComplexPolynomial& operator*=(double rhs);
        ComplexPolynomial& operator*=(const Complex& rhs);        
        ComplexPolynomial squared() const;
        Complex getCoefficient(unsigned int power) const;
        unsigned int getDegree() const;
        ComplexPolynomial getDerivative() const;
        bool operator==(const ComplexPolynomial& rhs) const;
        bool operator!=(const ComplexPolynomial& rhs) const;
        bool operator==(double rhs) const;
        bool operator!=(double rhs) const;
        bool operator==(int rhs) const;
        bool operator!=(int rhs) const;        
        friend std::ostream& operator<<(std::ostream& mycout, const ComplexPolynomial& p);

        // Slow distributive multiplication
        static ComplexPolynomial distributiveMultiplication(const ComplexPolynomial& lhs, const ComplexPolynomial& rhs);

        // Fast discrete fourier transform multiplication
        static ComplexPolynomial discreteFourierTransformMultiplication(ComplexPolynomial lhs, ComplexPolynomial rhs);

    protected:
        unsigned int degree_;
        std::vector<Complex> coefficients_;
        static double epsilon_; // Constant used for double comparison, doubles different by less than this amount are considered equal
        void unpad(); // Removes leading terms with a zero coefficient
        void pad(unsigned int degree); // Adds leading zeroes to make a polynomial seem like a specified degree
        bool static doubleComp(double a, double b);
        std::string static getPrintableComplex(const Complex& num);

        // /**
        //  * Evaluates the given polynomial on the nth roots of unity, recursively using FFT
        //  * 
        //  * @param f the function to evaluate
        //  * @param n which roots of unity to evaluate on, located in the roots of unity vector, must be a power of 2
        //  * @param rootsOfUnity a shared array of roots of unity to prevent them from being recalculated
        //  * @return a vector of the Complex numbers storing f evaluated at the roots of unity (in order going counter-clockwise around the unit circle)
        // */
        static std::vector<Complex> evaluateAtRootsOfUnity(const ComplexPolynomial& f, unsigned int n, const std::vector<std::vector<Complex>>& rootsOfUnity); 

        /**
         * Calculates the nth roots of unity
         * 
         * @param n the number of roots of unity to calculate
         * @return a vector of the roots of unity (ordered going counter-clockwise around the unit circle starting from 0)
        */
        static std::vector<Complex> calcRootsOfUnity(unsigned int n);

        static unsigned int getNextPowerOfTwo(unsigned int n);
};

/**
 * Adds an int to a Polynomial, the int will be treated as a degree 0 polynomial
 * 
 * @relatesalso Polynomial
*/
ComplexPolynomial operator+(int lhs, const ComplexPolynomial& rhs);

/**
 * Adds a double to a Polynomial, the double will be treated as a degree 0 polynomial
 * 
 * @relatesalso Polynomial
*/
ComplexPolynomial operator+(double lhs, const ComplexPolynomial& rhs);

ComplexPolynomial operator+(const Complex& lhs, const ComplexPolynomial& rhs);

/**
 * Subtracts an int to a Polynomial, the int will be treated as a degree 0 polynomial
 * 
 * @relatesalso Polynomial
*/
ComplexPolynomial operator-(int lhs, const ComplexPolynomial& rhs);

/**
 * Subtracts a double to a Polynomial, the double will be treated as a degree 0 polynomial
 * 
 * @relatesalso Polynomial
*/
ComplexPolynomial operator-(double lhs, const ComplexPolynomial& rhs);

ComplexPolynomial operator-(const Complex& lhs, const ComplexPolynomial& rhs);


/**
 * Multiplies a Polynomial by a constant
 * 
 * @relatesalso Polynomial
*/
ComplexPolynomial operator*(int lhs, const ComplexPolynomial& rhs);

/**
 * Multiplies a Polynomial by a constant
 * 
 * @relatesalso Polynomial
*/
ComplexPolynomial operator*(double lhs, const ComplexPolynomial& rhs);

ComplexPolynomial operator*(const Complex& lhs, const ComplexPolynomial& rhs);


bool operator==(double lhs, const ComplexPolynomial& rhs);
bool operator!=(double lhs, const ComplexPolynomial& rhs);
bool operator==(int lhs, const ComplexPolynomial& rhs);
bool operator!=(int lhs, const ComplexPolynomial& rhs);

template <typename V>
Complex ComplexPolynomial::evaluate(V x) const
{
    Complex answer = Complex();
    for(int i = degree_; i > 0; --i){ // View the polynomial in a factored form
        answer += coefficients_[i];
        answer *= x;
    }

    answer += coefficients_[0];

    return answer;
}

template <typename V>
Complex ComplexPolynomial::operator()(V x) const
{
    return evaluate(x);
}

#endif
