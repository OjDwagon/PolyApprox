#include "realPolynomial.h"
#include <vector>

using namespace std;

RealPolynomial::RealPolynomial():
    GenericPolynomial<double, RealPolynomial>()
{

}

RealPolynomial::RealPolynomial(const vector<double>& coefficients):
    GenericPolynomial<double, RealPolynomial>(coefficients)
{

}

RealPolynomial::RealPolynomial(const initializer_list<double>& coefficients):
    GenericPolynomial<double, RealPolynomial>(coefficients)
{

}

RealPolynomial RealPolynomial::operator*(const RealPolynomial& rhs) const
{
    return GenericPolynomial::distributiveMultiplication(*this, rhs);
}

RealPolynomial RealPolynomial::operator*(double rhs) const
{
    return *(static_cast<const GenericPolynomial<double, RealPolynomial>*>(this)) * rhs;
}

ostream& operator<<(ostream& out, const RealPolynomial& rhs)
{
    out << rhs.coefficients_[0];
    for(unsigned int i = 1; i <= rhs.degree_; ++i){
        out << " + " << rhs.coefficients_[i] << "x^" << i; 
    }
    return out;
}

RealPolynomial operator+(double lhs, const RealPolynomial& rhs)
{
    return rhs + lhs;
}

RealPolynomial operator-(double lhs, const RealPolynomial& rhs)
{
    return -1 * rhs + lhs;
}

RealPolynomial operator*(double lhs, const RealPolynomial& rhs)
{
    return static_cast<GenericPolynomial<double, RealPolynomial>>(rhs) * lhs;
}


bool operator==(double lhs, const RealPolynomial& rhs)
{
    return rhs == lhs;
}

bool operator!=(double lhs, const RealPolynomial& rhs)
{
    return rhs != lhs;
}
