#define _USE_MATH_DEFINES
#include <cmath>
#include "complex.h"

using namespace std;

double Complex::epsilon_ = 0.000001;

Complex::Complex(double primary, double other, bool cartesian):
    x_(), y_(), r_(), theta_(), cartesianStale_(!cartesian), polarStale_(cartesian)
{
    // Doesn't initialize the stale coords
    if(cartesian){
        x_ = primary;
        y_ = other;
    }
    else{
        r_ = primary;
        theta_ = other;
        normalizePolar();
    }
}

Complex::Complex(double xToBe, double yToBe, double rToBe, double thetaToBe):
    x_(xToBe), y_(yToBe), r_(rToBe), theta_(thetaToBe), cartesianStale_(false), polarStale_(false)
{
    normalizePolar();
}

Complex::Complex(double x):
    x_(x), y_(0), cartesianStale_(false), polarStale_(true)
{

}

double Complex::getCartesianX() const
{
    validateCartesian();    
    return x_;
}

double Complex::getCartesianY() const
{
    validateCartesian();
    return y_;
}

double Complex::getPolarR() const
{
    validatePolar();
    return r_;
}

double Complex::getPolarTheta() const
{
    validatePolar();
    return theta_;
}

void Complex::setCartesianX(double valToBe)
{
    validateCartesian(); // Validate Cartesian before trying to set
    x_ = valToBe;
    polarStale_ = true; // Don't update polar so it is now state
}

void Complex::setCartesianY(double valToBe)
{
    validateCartesian(); // Validate Cartesian before trying to set
    y_ = valToBe;
    polarStale_ = true; // Don't update polar so it is now state
}

void Complex::setPolarR(double valToBe)
{
    validatePolar(); // Validate Polar before trying to set
    r_ = valToBe;
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI
    cartesianStale_ = true; // Don't update cartesian so it is now state
}

void Complex::setPolarTheta(double valToBe)
{
    validatePolar(); // Validate Polar before trying to set
    theta_ = valToBe;
    normalizeTheta();
    cartesianStale_ = true; // Don't update cartesian so it is now state
}

Complex Complex::operator+(const Complex& rhs) const
{
    validateCartesian();
    rhs.validateCartesian();
    return Complex(x_ + rhs.x_, y_ + rhs.y_);
}

Complex Complex::operator+(double rhs) const
{
    validateCartesian();
    return Complex(x_ + rhs, y_);
}

Complex Complex::operator+(int rhs) const
{
    return *this + static_cast<double>(rhs);
}

Complex& Complex::operator+=(const Complex& rhs)
{
    // Validate Cartersian 
    validateCartesian();
    rhs.validateCartesian();

    // Perform arithmetic with cartesian
    x_ += rhs.x_;
    y_ += rhs.y_;

    // Polar is now stale
    polarStale_ = true;
    return *this;
}

Complex& Complex::operator+=(double rhs)
{
    // Validate Cartersian 
    validateCartesian();

    // Perform arithmetic with cartesian
    x_ += rhs;

    // Polar is now stale
    polarStale_ = true;
    return *this;
}

Complex& Complex::operator+=(int rhs)
{
    return *this += static_cast<double>(rhs);
}

Complex Complex::operator-(const Complex& rhs) const
{
    validateCartesian();
    rhs.validateCartesian();
    return Complex(x_ - rhs.x_, y_ - rhs.y_);
}

Complex Complex::operator-(double rhs) const
{
    validateCartesian();
    return Complex(x_ - rhs, y_);
}

Complex Complex::operator-(int rhs) const
{
    return *this - static_cast<double>(rhs);
}

Complex& Complex::operator-=(const Complex& rhs)
{
    // Validate Cartersian 
    validateCartesian();
    rhs.validateCartesian();

    // Perform arithmetic with cartesian
    x_ -= rhs.x_;
    y_ -= rhs.y_;

    // Polar is now stale
    polarStale_ = true;
    return *this;
}

Complex& Complex::operator-=(double rhs)
{
    // Validate Cartersian 
    validateCartesian();

    // Perform arithmetic with cartesian
    x_ -= rhs;

    // Polar is now stale
    polarStale_ = true;
    return *this;
}

Complex& Complex::operator-=(int rhs)
{
    return *this -= static_cast<double>(rhs);
}

Complex Complex::operator*(const Complex& rhs) const
{
    validatePolar();
    rhs.validatePolar();
    return Complex(r_ * rhs.r_, theta_ + rhs.theta_, false); // Do multiplication in polar form
}

Complex Complex::operator*(double rhs) const
{
    validateCartesian();
    validatePolar();
    return Complex(x_ * rhs, y_ * rhs, r_ * rhs, theta_); // Do multiplication in both cartesian and polar form
}

Complex Complex::operator*(int rhs) const
{
    return *this * static_cast<double>(rhs);
}

Complex& Complex::operator*=(const Complex& rhs)
{
    // Validate Polar coords
    validatePolar();
    rhs.validatePolar();

    // Perform operators in Polar coords
    r_ *= rhs.r_;
    theta_ += rhs.theta_;
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI

    // Cartesian is now stale
    cartesianStale_ = true;
    return *this;
}

Complex& Complex::operator*=(double rhs)
{
    // Perform operation in both coords for special case
    validateCartesian();
    validatePolar();

    x_ *= rhs;
    y_ *= rhs;
    r_ *= rhs;
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI
    return *this;
}

Complex& Complex::operator*=(int rhs)
{
    return *this *= static_cast<double>(rhs);
}

Complex Complex::operator/(const Complex& rhs) const
{
    validatePolar();
    rhs.validatePolar();
    return Complex(r_ / rhs.r_, theta_ - rhs.theta_, false); // Do division in polar form
}

Complex Complex::operator/(double rhs) const
{
    validateCartesian();
    validatePolar();

    return Complex(getCartesianX() / rhs, getCartesianY() / rhs, getPolarR() / rhs, getPolarTheta());
}

Complex Complex::operator/(int rhs) const
{
    return *this / static_cast<double>(rhs);
}

Complex& Complex::operator/=(const Complex& rhs)
{
    // Validate Polar coords
    validatePolar();
    rhs.validatePolar();

    // Perform operators in Polar coords
    r_ /= rhs.r_;
    theta_ -= rhs.theta_;
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI

    // Cartesian is now stale
    cartesianStale_ = true;
    return *this;
}

Complex& Complex::operator/=(double rhs)
{
    // Perform operation in both coords for special case
    validateCartesian();
    validatePolar();

    x_ /= rhs;
    y_ /= rhs;
    r_ /= rhs;
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI
    return *this;
}

Complex& Complex::operator/=(int rhs)
{
    return *this /= static_cast<double>(rhs);
}

Complex::operator double() const
{
    validateCartesian();
    return x_;
}

const Complex& Complex::freshen() const
{
    validateCartesian();
    validatePolar();
    return *this;
}

void Complex::validateCartesian() const
{
    if(cartesianStale_){
        updateCartesian();
    }
}

void Complex::validatePolar() const
{
    if(polarStale_){
        updatePolar();
    }
}

void Complex::updateCartesian() const
{
    x_ = calcCartesianX();
    y_ = calcCartesianY();
    cartesianStale_ = false;
}

void Complex::updatePolar() const
{
    r_ = calcPolarR();
    theta_ = calcPolarTheta();
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI;
    polarStale_ = false;
}

double Complex::calcCartesianX() const
{
    return r_ * cos(theta_);
}

double Complex::calcCartesianY() const
{
    return r_ * sin(theta_);
}

double Complex::calcPolarR() const
{
    return sqrt(x_ * x_ + y_ * y_);
}

double Complex::calcPolarTheta() const
{
    if(doubleComp(x_, 0)){ // arctan is undefined when x is 0
        if(doubleComp(y_, 0)){ // Special case of (0, 0)
            return 0;
        }
        else if(y_ > 0){ // Aligned with positive y axis
            return M_PI / 2;
        }
        else{ // Aligned with negative y axis
            return 1.5 * M_PI;
        }
    }
    else{
        // Calculate theta by quadrant
        double theta = atan2(y_, x_); // Gives a range from -PI to PI
        if(theta < 0){ // Converts to a range from 0 to 2PI
            theta += 2 * M_PI;
        }
        return theta;
    }
}

void Complex::normalizePolar() const
{
    if(r_ < 0){
        r_ *= -1;
        theta_ += M_PI;
    }
    normalizeTheta();
}

void Complex::normalizeTheta() const
{
    theta_ = fmod(theta_, 2 * M_PI);
}

Complex operator+(double lhs, const Complex& rhs)
{
    return rhs + lhs;
}

Complex operator+(int lhs, const Complex& rhs)
{
    return rhs + lhs;
}

Complex operator-(double lhs, const Complex& rhs)
{
    return -1 * rhs + lhs;
}

Complex operator-(int lhs, const Complex& rhs)
{
    return -1 * rhs + lhs;
}

Complex operator*(double lhs, const Complex& rhs)
{
    return rhs * lhs;
}

Complex operator*(int lhs, const Complex& rhs)
{
    return rhs * lhs;
}

Complex operator/(double lhs, const Complex& rhs)
{
    return Complex(lhs, 0) / rhs;
}

Complex operator/(int lhs, const Complex& rhs)
{
    return Complex(lhs, 0) / rhs;
}

Complex Complex::operator^(double n) const
{
    validatePolar();
    return Complex(pow(r_, n), theta_ * n, false);
}

Complex Complex::operator^(int n) const
{
    return *this ^ static_cast<double>(n);
}

Complex& Complex::operator^=(double n)
{
    validatePolar();
    r_ = pow(r_, n);
    theta_ *= n;
    normalizePolar(); // Normalize polar so r is positive and theta is from 0 to 2PI

    cartesianStale_ = true;
    return *this;
}

Complex& Complex::operator^=(int n)
{
    return *this ^= static_cast<double>(n);
}

bool Complex::operator==(const Complex& rhs) const
{
    validateCartesian();
    rhs.validateCartesian();

    return doubleComp(x_, rhs.x_) && doubleComp(y_, rhs.y_); 
}

bool Complex::operator!=(const Complex& rhs) const
{
    return !(*this == rhs);
}

bool Complex::operator==(double rhs) const
{
    validateCartesian();
    
    return doubleComp(x_, rhs) && doubleComp(y_, 0);
}

bool Complex::operator!=(double rhs) const
{
    return !(*this == rhs);
}

bool Complex::operator==(int rhs) const
{
    return *this == static_cast<double>(rhs); // Just cast to double and do comparison that way
}

bool Complex::operator!=(int rhs) const
{
    return !(*this == rhs);
}

bool operator==(double lhs, const Complex& rhs)
{
    return rhs == lhs;
}

bool operator!=(double lhs, const Complex& rhs)
{
    return !(rhs == lhs);
}

bool operator==(int lhs, const Complex& rhs)
{
    return rhs == lhs;
}

bool operator!=(int lhs, const Complex& rhs)
{
    return !(rhs == lhs);
}

bool Complex::isReal() const
{
    return doubleComp(y_, 0);
}

bool Complex::doubleComp(double a, double b)
{
    return ((a - b) < epsilon_) && ((b - a) < epsilon_);
}

void Complex::charactarization(double& outX, double& outY, double& outR, double& outTheta, bool& outCartesianStale, bool& outPolarStale) const
{
    outX = x_;
    outY = y_;
    outR = r_;
    outTheta = theta_;
    outCartesianStale = cartesianStale_;
    outPolarStale = polarStale_;
}

ostream& operator<<(ostream& mycout, const Complex& rhs)
{
    if(rhs.cartesianStale_){
        mycout << "(stale, ";
    }
    else{
        mycout << "(fresh,  ";
    }
    mycout << rhs.x_ << " +  " << rhs.y_ << "i) ";

    if(rhs.polarStale_){
        mycout << "(stale, ";
    }
    else{
        mycout << "(fresh, ";
    }
    mycout << "r= " << rhs.r_ << ", theta = " << rhs.theta_ << ")";

    return mycout;
}


