#include "realPolynomial.h"
#include "complexPolynomial.h"
#include <vector>
#include <functional>

using namespace std;

double RealPolynomial::epsilon_ = 0.000001;
function <bool(double, double)> RealPolynomial::comp = [](double a, double b) -> bool {return ((a - b) < epsilon_) && ((b - a) < epsilon_);};

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
    unsigned int minDegree = min(degree_, rhs.degree_);
    unsigned int diff = rhs.degree_ + degree_ -  2 * minDegree;
    if(minDegree >= 4096 && diff < minDegree / 2 ){ // Use fourier transform division when the polynomials are balanced and of a high degree
        // Convert to complex polynomials, perform multiplication, then return answer as real polynomial
        return discreteFourierTransformMultiplication(*this, rhs);
    }
    else{ // Else use distributive multiplication
        return distributiveMultiplication(*this, rhs);
    }
}

RealPolynomial RealPolynomial::operator*(double rhs) const
{
    return *(static_cast<const GenericPolynomial<double, RealPolynomial>*>(this)) * rhs;
}

RealPolynomial RealPolynomial::discreteFourierTransformMultiplication(const RealPolynomial& lhs, const RealPolynomial& rhs)
{
    // Convert to complex polynomials, perform multiplication, then return answer as real polynomial
    ComplexPolynomial complexLHS = ComplexPolynomial(lhs.getCoefficients());
    ComplexPolynomial complexRHS = ComplexPolynomial(rhs.getCoefficients());
    vector<Complex> product = ComplexPolynomial::discreteFourierTransformMultiplication(complexLHS, complexRHS).getCoefficients();
    return RealPolynomial(vector<double>(product.begin(), product.end()));
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
