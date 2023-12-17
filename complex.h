#ifndef COMPLEX_H
#define COMPLEX_H

#include <iostream>

/**
 * A complex number class that models complex numbers and their operations,
 * storing as both cartesian and polar representation,
 * Note: instead of an i component, a y component is used instead
 * 
 * A lazy implementation that only updates cartesian or polar form rendering the other stale
*/
class Complex
{
    public:
        /**
         * Constructor for a Complex number given only 1 type of coordinates, defaults to 0
         * The other type of coordiantes is left uninitialized and is stale
         * 
         * @param r the real(x) component or radius component (depending on if in cartesian or polar)
         * @param other the imaginary(y) or angle componenent (depending on if in cartesian or polar)
         * @param cartesian a bool to indicate if coordinates are provided in cartesain or polar form
        */
        Complex(double primary = 0, double other = 0, bool cartesian = true);

        /**
         * Constructor for a complex number given 2 types of coordiantes
         * 
         * @param xToBe the x component in cartesian form
         * @param yToBe the y component in cartesian form
         * @param rToBe the r value in polar form
         * @param thetaToBe the theta value in polar form
        */
        Complex(double xToBe, double yToBe, double rToBe, double thetaToBe);

        /**
         * Getter for x cartesian coordinate, will update cartesian coordinates if stale for faster operations later
        */
        double getCartesianX() const;

        /**
         * Getter for y cartesian coordinate, will update cartesian coordinates if stale for faster operations later
        */
        double getCartesianY() const;

        /**
         * Getter for r polar coordinate, will update polar coordinates if stale for faster operations later
        */
        double getPolarR() const;

        /**
         * Getter for the theta polar coordinate, will update polar coordinates if stale for faster operations later
        */
        double getPolarTheta() const;

        /**
         * Setter for the x cartesian coordinate, makes polar coords stale
        */
        void setCartesianX(double valToBe);

        /**
         * Setter for the y cartesian coordinate, makes polar coords stale
        */
        void setCartesianY(double valToBe);

        /**
         * Setter for the r polar coordinate, makes cartesian coords stale
        */
        void setPolarR(double valToBe);

        /**
         * Setter for the theta polar coordinate, makes cartesian coords stale
        */
        void setPolarTheta(double valToBe);

        /**
         * Adds 2 Complex numbers using cartesian coodinates can interpret doubles as complex numbers
         * Does not update polar coordinates making them stale
        */
        Complex operator+(const Complex& rhs) const;
        Complex operator+(double rhs) const;
        Complex operator+(int rhs) const;
        Complex& operator+=(const Complex& rhs);
        Complex& operator+=(double rhs);
        Complex& operator+=(int rhs);

        /**
         * Subtracts 2 Complex numbers using cartesian coodinates
         * Does not update polar coordinates making them stale
        */
        Complex operator-(const Complex& rhs) const;
        Complex operator-(double rhs) const;
        Complex operator-(int rhs) const;
        Complex& operator-=(const Complex& rhs);
        Complex& operator-=(double rhs);
        Complex& operator-=(int rhs);

        /**
         * Multiplies 2 Complex numbers using polar coordinates
         * Does not update cartesian coordinates making them stale
        */
        Complex operator*(const Complex& rhs) const;

        /**
         * Multiplies a Complex number by a constant
         * Updates both cartesian and polar coordinates
        */
        Complex operator*(double rhs) const;
        Complex operator*(int rhs) const;
        Complex& operator*=(const Complex& rhs);
        Complex& operator*=(double rhs);
        Complex& operator*=(int rhs);

        /**
         * Divides 2 Complex numbers using polar coordinates
         * Does not update cartesian coordinates making them stale
        */
        Complex operator/(const Complex& rhs) const;
        Complex operator/(double rhs) const;
        Complex operator/(int rhs) const;
        Complex& operator/=(const Complex& rhs);
        Complex& operator/=(double rhs);
        Complex& operator/=(int rhs);

        /**
         * Raises a complex number to a power
         * 
         * @param n a real value number to raise to the power of
        */
        Complex operator^(double n) const;
        Complex operator^(int n) const;
        Complex& operator^=(double n);
        Complex& operator^=(int n);

        /**
         * Performs comparisons in cartesian coords
        */
        bool operator==(const Complex& rhs) const;
        bool operator!=(const Complex& rhs) const;
        bool operator==(double rhs) const;
        bool operator!=(double rhs) const;
        bool operator==(int rhs) const;
        bool operator!=(int rhs) const;

        operator double() const;
        
        /**
         * Refreshes cartesian or polar coordinates so that neither are stale and returns a reference to it
        */
        const Complex& freshen() const;

        bool isReal() const;

        friend std::ostream& operator<<(std::ostream& mycout, const Complex& rhs);

        /**
         * Debug function that returns all the numbers that characterize a complex number and its state
        */
        void charactarization(double& outX, double& outY, double& outR, double& outTheta, bool& outCartesianStale, bool& outPolarStale) const;
    protected:
        mutable double x_;
        mutable double y_;
        mutable double r_; // A positive (or zero) value
        mutable double theta_; // In radians between [0, 2PI)
        mutable bool cartesianStale_;
        mutable bool polarStale_;
        static double epsilon_;


        void validateCartesian() const;
        void validatePolar() const;
        void updateCartesian() const;
        void updatePolar() const;
        double calcCartesianX() const;
        double calcCartesianY() const;
        double calcPolarR() const;
        double calcPolarTheta() const;
        void normalizePolar() const; // Normalize r then theta in polar form
        void normalizeTheta() const; 

        bool doubleComp(double a, double b) const;
   

};

Complex operator+(double lhs, const Complex& rhs);
Complex operator+(int lhs, const Complex& rhs);
Complex operator-(double lhs, const Complex& rhs);
Complex operator-(int lhs, const Complex& rhs);
Complex operator*(double lhs, const Complex& rhs);
Complex operator*(int lhs, const Complex& rhs);
Complex operator/(double lhs, const Complex& rhs);
Complex operator/(int lhs, const Complex& rhs);
bool operator==(double lhs, const Complex& rhs);
bool operator!=(double lhs, const Complex& rhs);
bool operator==(int lhs, const Complex& rhs);
bool operator!=(int lhs, const Complex& rhs);


/**
 * Outputs a representation of a Complex number to an output stream
 * 
 * @relatesalso Complex
*/
std::ostream& operator<<(std::ostream& mycout, const Complex& rhs);


#endif
