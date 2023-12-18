
#include "complexPolynomial.h"
#include "complex.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

double ComplexPolynomial::epsilon_ = 0.000001;

ComplexPolynomial::ComplexPolynomial():
    degree_(0), coefficients_({0})
{

}

ComplexPolynomial::ComplexPolynomial(const vector<Complex>& coefficients):
    degree_(coefficients.size() - 1), coefficients_(coefficients)
{
    unpad();

}

void ComplexPolynomial::setCoefficient(unsigned int power, Complex val)
{
    if(power <= degree_){ // If it is a valid coefficient
        coefficients_[power] = val;
    }
    else{ // If it needs to be expanded 
        while(coefficients_.size() < power){
            coefficients_.push_back(Complex());
        }
        coefficients_.push_back(val);
        degree_ = power;
    }

    unpad();
}

ComplexPolynomial ComplexPolynomial::operator+(const ComplexPolynomial& rhs) const
{
    // Assumes lhs is bigger, if rhs is bigger flip them
    if(rhs.degree_ > degree_){
        return rhs + *this;
    }

    ComplexPolynomial lhs = *this; // Make a copy

    for(unsigned int i = 0; i <= rhs.degree_; ++i){ // Add all the coefficients
        lhs.coefficients_[i] += rhs.coefficients_[i];
    }

    return lhs; // Return it
}

ComplexPolynomial ComplexPolynomial::operator+(int rhs) const
{
    ComplexPolynomial lhs = *this; // Make a copy
    lhs.coefficients_[0] += rhs; // Add the constant
    return lhs; // Return it
}

ComplexPolynomial ComplexPolynomial::operator+(double rhs) const
{
    ComplexPolynomial lhs = *this; // Make a copy
    lhs.coefficients_[0] += rhs; // Add the constant
    return lhs; // Return it
}

ComplexPolynomial ComplexPolynomial::operator+(const Complex& rhs) const
{
    ComplexPolynomial lhs = *this; // Make a copy
    lhs.coefficients_[0] += rhs; // Add the constant
    return lhs; // Return it
}

ComplexPolynomial& ComplexPolynomial::operator+=(const ComplexPolynomial& rhs)
{
    if(rhs.degree_ > degree_){ // Pad up if necessary
        pad(rhs.degree_);
    }

    for(unsigned int i = 0; i <= rhs.degree_; ++i){ // Add all the coefficients
        coefficients_[i] += rhs.coefficients_[i];
    }

    unpad();

    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator+=(int rhs)
{
    coefficients_[0] += rhs; // Add the constant
    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator+=(double rhs)
{
    coefficients_[0] += rhs; // Add the constant
    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator+=(const Complex& rhs)
{
    coefficients_[0] += rhs; // Add the constant
    return *this;
}

ComplexPolynomial ComplexPolynomial::operator-(const ComplexPolynomial& rhs) const
{
    return *this + -1 * rhs;
}

ComplexPolynomial ComplexPolynomial::operator-(int rhs) const
{
    return *this + -1 * rhs;
}

ComplexPolynomial ComplexPolynomial::operator-(double rhs) const
{
    return *this + -1 * rhs;
}

ComplexPolynomial ComplexPolynomial::operator-(const Complex& rhs) const
{
    return *this + -1 * rhs;
}

ComplexPolynomial& ComplexPolynomial::operator-=(const ComplexPolynomial& rhs)
{
    *this += -1 * rhs;
    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator-=(int rhs)
{
    coefficients_[0] -= rhs; // Subtract the constant
    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator-=(double rhs)
{
    coefficients_[0] -= rhs; // Subtract the constant
    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator-=(const Complex& rhs)
{
    coefficients_[0] -= rhs; // Subtract the constant
    return *this;
}

ComplexPolynomial ComplexPolynomial::operator*(const ComplexPolynomial& rhs) const
{
    return ComplexPolynomial(); // Placeholder
}

ComplexPolynomial ComplexPolynomial::operator*(int rhs) const
{
    vector<Complex> product = this->coefficients_; // Make a copy

    for(unsigned int i = 0; i <= this->degree_; ++i){ // Multiply everything
        product[i] *= rhs;
    }

    return ComplexPolynomial(product); // Return
}

ComplexPolynomial ComplexPolynomial::operator*(double rhs) const
{
    vector<Complex> product = this->coefficients_; // Make a copy

    for(unsigned int i = 0; i <= this->degree_; ++i){ // Multiply everything
        product[i] *= rhs;
    }


    return ComplexPolynomial(product); // Return
}

ComplexPolynomial ComplexPolynomial::operator*(const Complex& rhs) const
{
    vector<Complex> product = this->coefficients_; // Make a copy

    for(unsigned int i = 0; i <= this->degree_; ++i){ // Multiply everything
        product[i] *= rhs;
    }


    return ComplexPolynomial(product); // Return
}

ComplexPolynomial& ComplexPolynomial::operator*=(const ComplexPolynomial& rhs)
{
    return *this; // Placeholder
}

ComplexPolynomial& ComplexPolynomial::operator*=(int rhs)
{
    for(unsigned int i = 0; i <= degree_; ++i){
        coefficients_[i] *= rhs;
    }

    unpad();

    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator*=(double rhs)
{
    for(unsigned int i = 0; i <= degree_; ++i){
        coefficients_[i] *= rhs;
    }

    unpad();

    return *this;
}

ComplexPolynomial& ComplexPolynomial::operator*=(const Complex& rhs)
{
    for(unsigned int i = 0; i <= degree_; ++i){
        coefficients_[i] *= rhs;
    }

    unpad();

    return *this;
}


bool ComplexPolynomial::operator==(const ComplexPolynomial& rhs) const
{
    // Two polynomials of different degree can't be equal
    if(degree_ != rhs.degree_){
        return false;
    }

    // Check if all coefficients match
    for(unsigned int i = 0; i < coefficients_.size(); ++i){
        if(coefficients_[i] != rhs.coefficients_[i]){
            return false;
        }
    }

    return true;
    
}

bool ComplexPolynomial::operator!=(const ComplexPolynomial& rhs) const
{
    return !(*this == rhs);
}

bool ComplexPolynomial::operator==(double rhs) const
{
    if(degree_ == 0){
        return rhs == coefficients_[0];
    }

    return false;
}

bool ComplexPolynomial::operator!=(double rhs) const
{
    return !(*this == rhs);
}

bool ComplexPolynomial::operator==(int rhs) const
{
    if(degree_ == 0){
        return rhs == coefficients_[0];
    }

    return false;
}

bool ComplexPolynomial::operator!=(int rhs) const
{
    return !(*this == rhs);
}

bool operator==(double lhs, const ComplexPolynomial& rhs)
{
    return rhs == lhs;
}

bool operator!=(double lhs, const ComplexPolynomial& rhs)
{
    return rhs != lhs;
}

bool operator==(int lhs, const ComplexPolynomial& rhs)
{
    return rhs == lhs;
}

bool operator!=(int lhs, const ComplexPolynomial& rhs)
{
    return rhs != lhs;
}



ComplexPolynomial ComplexPolynomial::squared() const
{
    return (*this) * (*this);
}

Complex ComplexPolynomial::getCoefficient(unsigned int power) const
{
    if(power <= degree_){
        return coefficients_[power];
    }
    else{
        return 0;
    }
}

unsigned int ComplexPolynomial::getDegree() const
{
    return degree_;
}

ComplexPolynomial ComplexPolynomial::getDerivative() const
{
    vector<Complex> derivative(degree_);

    for(unsigned int i = 0; i < degree_; ++i){
        derivative[i] = (int)(i + 1) * coefficients_[i + 1];
    }

    return ComplexPolynomial(derivative);
}

ComplexPolynomial ComplexPolynomial::distributiveMultiplication(const ComplexPolynomial& lhs, const ComplexPolynomial& rhs)
{
    ComplexPolynomial product;
    product.pad(lhs.degree_ + rhs.degree_);

    for(unsigned int i = 0; i <= lhs.degree_; ++i){
        for(unsigned int j = 0; j <= rhs.degree_; ++j){
            product.coefficients_[i + j] += lhs.coefficients_[i] * rhs.coefficients_[j];
        }
    }

    return product;
}

static ComplexPolynomial discreteFourierTransformMultiplication(const ComplexPolynomial& lhs, const ComplexPolynomial& rhs)
{
    
}

void ComplexPolynomial::unpad()
{
    int counter = 0;
    for(unsigned int i = degree_; i > 0; --i){
        if(coefficients_[i] == Complex()){
            ++counter;
        }
        else{
            break;
        }
    }

    degree_ -= counter;
    for(int i = 0; i < counter; ++i){
        coefficients_.pop_back();
    }
}

void ComplexPolynomial::pad(unsigned int degree)
{
    if(degree_ >= degree){ // Do nothing if you are already at or above the required degree
        return;
    }

    for(unsigned int i = 0; i < degree - degree_; ++i){ // Add zeroes to the back
        coefficients_.push_back(Complex());
    }

    degree_ = degree; // Update the degree
}

bool ComplexPolynomial::doubleComp(double a, double b)
{
    return ((a - b) < epsilon_) && ((b - a) < epsilon_);
}

string ComplexPolynomial::getPrintableComplex(const Complex& num)
{
    stringstream ss;
    ss << "(" << num.getCartesianX() << " + " << num.getCartesianY() << "i)";
    return ss.str();
}

ComplexPolynomial operator+(int lhs, const ComplexPolynomial& rhs)
{
    return rhs + lhs;
}

ComplexPolynomial operator+(double lhs, const ComplexPolynomial& rhs)
{
    return rhs + lhs;
}

ComplexPolynomial operator+(const Complex& lhs, const ComplexPolynomial& rhs)
{
    return rhs + lhs;
}

ComplexPolynomial operator-(int lhs, const ComplexPolynomial& rhs)
{
    return -1 * rhs + lhs;
}

ComplexPolynomial operator-(double lhs, const ComplexPolynomial& rhs)
{
    return -1 * rhs + lhs;
}

ComplexPolynomial operator-(const Complex& lhs, const ComplexPolynomial& rhs)
{
    return -1 * rhs + lhs;
}

ComplexPolynomial operator*(int lhs, const ComplexPolynomial& rhs)
{
    return rhs * lhs;
}

ComplexPolynomial operator*(double lhs, const ComplexPolynomial& rhs)
{
    return rhs * lhs;
}

ComplexPolynomial operator*(const Complex& lhs, const ComplexPolynomial& rhs)
{
    return rhs * lhs;
}

ostream& operator<<(ostream& mycout, const ComplexPolynomial& p)
{
    mycout << ComplexPolynomial::getPrintableComplex(p.coefficients_[0]);
    
    for(unsigned int i = 1; i <= p.degree_; ++i){
        mycout << " + " << p.getPrintableComplex(p.coefficients_[i]) << "x^" << i;
    }

    return mycout;
}



