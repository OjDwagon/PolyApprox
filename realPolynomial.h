#ifndef REAL_POLYNOMIAL_H
#define REAL_POLYNOMIAL_H

#include <vector>
#include "genericPolynomial.h"

#include <iostream>

/**
 * A class the models real valued polynomials of a non-negative degree and their operations (excluding division)
*/
class RealPolynomial: public GenericPolynomial<double, RealPolynomial>
{
    public:
        RealPolynomial();
        RealPolynomial(const std::vector<double>& coefficients);
        RealPolynomial(const std::initializer_list<double>& coefficients);
        RealPolynomial operator*(const RealPolynomial& rhs) const;
        RealPolynomial operator*(double rhs) const;
        friend std::ostream& operator<<(std::ostream& out, const RealPolynomial& rhs);

    protected:
    
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
