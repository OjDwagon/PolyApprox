#include "realPolynomial.h"
#include "complexPolynomial.h"
#include <vector>

using namespace std;

double RealPolynomial::epsilon_ = 0.000001;

RealPolynomial::RealPolynomial():
    GenericPolynomial<double, RealPolynomial>()
{
    // Set the function's comparator to use double comparison instead of the default equality
    GenericPolynomial<double, RealPolynomial>::comp = [](double a, double b) -> bool {return ((a - b) < epsilon_) && ((b - a) < epsilon_);};
}

RealPolynomial::RealPolynomial(const vector<double>& coefficients):
    GenericPolynomial<double, RealPolynomial>(coefficients)
{
    // Set the function's comparator to use double comparison instead of the default equality
    GenericPolynomial<double, RealPolynomial>::comp = [](double a, double b) -> bool {return ((a - b) < epsilon_) && ((b - a) < epsilon_);};
}

RealPolynomial::RealPolynomial(const initializer_list<double>& coefficients):
    GenericPolynomial<double, RealPolynomial>(coefficients)
{
    // Set the function's comparator to use double comparison instead of the default equality
    GenericPolynomial<double, RealPolynomial>::comp = [](double a, double b) -> bool {return ((a - b) < epsilon_) && ((b - a) < epsilon_);};
}

RealPolynomial RealPolynomial::operator*(const RealPolynomial& rhs) const
{
    return GenericPolynomial::distributiveMultiplication(*this, rhs);
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
