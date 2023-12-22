#ifndef COMPLEX_POLYNOMIAL_H
#define COMPLEX_POLYNOMIAL_H

#include <vector>
#include "complex.h"
#include <string>
#include <iostream>

/**
 * A class the models complex valued polynomials of a non-negative degree and their operators (excluding division)
*/
class ComplexPolynomial 
{
    public:
        /**
         * Default constructor for a ComplexPolynomial, defaults to the constant 0
        */
        ComplexPolynomial();

        /**
         * ComplexPolynomial constructor given coefficients
         * 
         * @tparam V the type of the coefficients, must support casting to a Complex number
         * @param coefficients a vector of coefficeints, the ith element of the vector is used as the coefficient for x^i
        */
        template <typename V>
        ComplexPolynomial(const std::vector<V>& coefficients);

        /**
         * ComplexPolynomial constructor given coefficients
         * 
         * @tparam V the type of the coefficients, must support casting to a Complex number
         * @param coefficients an initinialization of complex numbers, the ith element of the vector is used as the coefficient for x^i
        */
        template <typename V>
        ComplexPolynomial(const std::initializer_list<V>& coefficients);       

        /**
         * Sets the coefficient corresponding to x^power to be val, does nothing if a negative power is given
         * 
         * @param power the degree of x corresponding to the coefficient to be set
         * @param val the value that the coefficient will be set to
        */
        void setCoefficient(unsigned int power, Complex val); // Does nothing if an invalid power is given

        /**
         * A linear time evaluation method of the polynomial
         * 
         * @tparam V the type of the input into the polynomial, must have defined arithmetic operations with Complex Numbers
         * @param x the value that the polynomial will be evaluated at
         * 
         * @returns the value of the polynomial evaluated at x as a complex number
        */
        template <typename V>
        Complex evaluate(V x) const;

        /**
         * Operator() that allows for evaluation of a polynomial in function style syntax
         * 
         * @tparam V the type of the input into the polynomial, must have defined arithmetic operations with Complex Numbers
         * @param x the value that the polynomial will be evaluated at
         * 
         * @returns the value of the polynomial evaluated at x as a complex number
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

        /**
         * Calculates and returns the square of this polynomial
         * 
         * @returns a ComplexPolynomial which is the result of squaring this polynomial
        */
        ComplexPolynomial squared() const;

        /**
         * Calculates and returns the polynommial raised to the desired non-negative power
         * 
         * @param power a non-negative number which is the power the polynomial will be raised to
         * 
         * @returns a ComplexPolynomial which is the result of raising this polynomial to the specified power
        */
        ComplexPolynomial exponentiate(unsigned int power) const;

        /**
         * Getter for the coefficient associated with x^power, returns 0 if given an invalid coefficient
         * 
         * @param power the degree of x corresponding to the coefficient to be set
         * 
         * @returns the coefficient associated with x^power
        */
        Complex getCoefficient(unsigned int power) const;

        /**
         * Getter for the degree of the polynomial, which is the highest n such that the coefficient for x^n is not 0
         * 
         * @returns the degree of the ComplexPolynomial
        */
        unsigned int getDegree() const;

        /**
         * Applies the power rule to calculate the derivative of this polynomial
         * 
         * @returns a ComplexPolynomial which is the derivate of this polynomial
        */
        ComplexPolynomial getDerivative() const;


        bool operator==(const ComplexPolynomial& rhs) const;
        bool operator!=(const ComplexPolynomial& rhs) const;
        bool operator==(double rhs) const;
        bool operator!=(double rhs) const;
        bool operator==(int rhs) const;
        bool operator!=(int rhs) const;        
        friend std::ostream& operator<<(std::ostream& mycout, const ComplexPolynomial& p);

        /**
         * Returns the product of 2 polynomials, implemented using standard distribution in quadratic time
         * 
         * @param lhs the first polynomial to multiply
         * @param rhs the second polynomial to multiply
         * 
         * @returns the polynomial lhs * rhs
        */
        static ComplexPolynomial distributiveMultiplication(const ComplexPolynomial& lhs, const ComplexPolynomial& rhs);

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
        unsigned int degree_;
        std::vector<Complex> coefficients_;
        static double epsilon_; // Constant used for double comparison, doubles different by less than this amount are considered equal

        /**
         * Removes all leading zero coefficients (ignoring the constant coefficent x^0)
        */
        void unpad();

        /**
         * Adds leading zero coefficients to a polynomial to make it seem like a specified degree
         * 
         * @param degree the desired degree of the polynomial after the padding operation
        */
        void pad(unsigned int degree);

        /**
         * Helper function for double comparison, returns true if the doubels are within epsilon_ of each other
         * 
         * @param a, b the numbers to be compared
         * 
         * @returns true if a and b are within epsilon_ of each other
        */
        bool static doubleComp(double a, double b);

        /**
         * Returns a printable string of a Complex number in standard form (cartesian coordinates)
         * 
         * @param num the Complex number to be represented
         * @returns a string representing the Complex number
        */
        std::string static getPrintableComplex(const Complex& num);

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
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator+(double lhs, const ComplexPolynomial& rhs);

/**
 * Adds a Complex number to a Polynomial, the Complex number will be treated as a degree 0 polynomial
 * 
 * @relatesalso ComplexPolynomial
*/
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
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator-(double lhs, const ComplexPolynomial& rhs);

/**
 * Subtracts a Complex number to a Polynomial, the Complex number will be treated as a degree 0 polynomial
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator-(const Complex& lhs, const ComplexPolynomial& rhs);


/**
 * Multiplies a Polynomial by a constant
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator*(int lhs, const ComplexPolynomial& rhs);

/**
 * Multiplies a Polynomial by a constant
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator*(double lhs, const ComplexPolynomial& rhs);


/**
 * Multiplies a Polynomial by a constant
 * 
 * @relatesalso ComplexPolynomial
*/
ComplexPolynomial operator*(const Complex& lhs, const ComplexPolynomial& rhs);


bool operator==(double lhs, const ComplexPolynomial& rhs);
bool operator!=(double lhs, const ComplexPolynomial& rhs);
bool operator==(int lhs, const ComplexPolynomial& rhs);
bool operator!=(int lhs, const ComplexPolynomial& rhs);

template <typename V>
ComplexPolynomial::ComplexPolynomial(const std::vector<V>& coefficients):
    degree_(coefficients.size() - 1), coefficients_(coefficients.begin(), coefficients.end())
{
    unpad();
}

template <typename V>
ComplexPolynomial::ComplexPolynomial(const std::initializer_list<V>& coefficients):
    degree_(coefficients.size() - 1)
{
    for(const V& c: coefficients){
        coefficients_.push_back(Complex(c));
    }

    unpad();
}

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
