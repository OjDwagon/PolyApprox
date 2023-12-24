#ifndef REAL_POLYNOMIAL_H
#define REAL_POLYNOMIAL_H

#include <vector>
#include "genericPolynomial.h"

#include <iostream>
#include <functional>

/**
 * @brief A class the models real valued polynomials of a non-negative degree and their operations (excluding division)
*/
class RealPolynomial: public GenericPolynomial<double, RealPolynomial>
{
    public:
        /**
         * Default constructor for a RealPolynomial, defaults to constant 0
        */
        RealPolynomial();

        /**
         * RealPolynomial constructor given coefficients
         * 
         * @param coefficients a vector of coefficeints, the ith element of the vector is used as the coefficient for x^i
        */
        RealPolynomial(const std::vector<double>& coefficients);

        /**
         * RealPolynomial constructor given coefficients
         * 
         * @param coefficients an initialization list of type double, the ith element of the list is used as the coefficient for x^i
        */
        RealPolynomial(const std::initializer_list<double>& coefficients);
        
        RealPolynomial operator*(const RealPolynomial& rhs) const;
        RealPolynomial operator*(double rhs) const;
        friend std::ostream& operator<<(std::ostream& out, const RealPolynomial& rhs);

        /**
         * Returns the product of 2 polynomials, implemented using the Discrete Fourier Transform via ComplexPolynomials in loglinear time
         * 
         * @param lhs the first polynomial to multiply
         * @param rhs the second polynomial to multiply
         * 
         * @returns the polynomial lhs * rhs
        */
        static RealPolynomial discreteFourierTransformMultiplication(const RealPolynomial& lhs, const RealPolynomial& rhs);

    protected:
        static double epsilon_; // Constant used for double comparison, doubles different by less than this amount are considered equal
        friend GenericPolynomial<double, RealPolynomial>;
        static std::function<bool(double, double)> comp;
};

/**
 * Adds a double to a RealPolynomial, the double will be treated as a degree 0 RealPolynomial
 * 
 * @relatesalso RealPolynomial
*/
RealPolynomial operator+(double lhs, const RealPolynomial& rhs);

/**
 * Subtracts a double to a RealPolynomial, the double will be treated as a degree 0 RealPolynomial
 * 
 * @relatesalso RealPolynomial
*/
RealPolynomial operator-(double lhs, const RealPolynomial& rhs);

/**
 * Multiplies a RealPolynomial by a constant
 * 
 * @relatesalso RealPolynomial
*/
RealPolynomial operator*(double lhs, const RealPolynomial& rhs);

/**
 * Equality check of a RealPolynomial with a constant, returns true if the RealPolynomial is degree 0 and has a matching coefficient
 * 
 * @relatesalso RealPolynomial
*/
bool operator==(double lhs, const RealPolynomial& rhs);

/**
 * Inequality check of a RealPolynomial with a constant, returns false if the RealPolynomial is degree 0 and has a matching coefficient
 * 
 * @relatesalso RealPolynomial
*/
bool operator!=(double lhs, const RealPolynomial& rhs);


#endif
