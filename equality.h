#ifndef EQUALITY_H
#define EQUALITY_H
#include <cmath>
#include "complex.h"

class Complex;

/**
 * A functor for double comparison
*/
struct DoubleEquality
{
    public:
        DoubleEquality(double epsilon = 0.000001);

        /**
         * Returns true if two doubles are within an error of each other
         * 
         * @param x the first double
         * @param y the second double
        */
        bool operator()(double x, double y);
        bool operator()(const Complex& x, const Complex& y);

    private:
        double epsilon_; // The error allowed between two doubles considered the same
};


template <typename T>
struct StandardEquality
{
    bool operator()(T x, T y)
    {
        return x == y;
    }
};

#endif
