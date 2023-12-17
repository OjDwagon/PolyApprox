
#include "equality.h"
#include "complex.h"
#include <iostream>

DoubleEquality::DoubleEquality(double epsilon):
    epsilon_(epsilon)
{

}

bool DoubleEquality::operator()(double x, double y)
{
    return ((x - y) < epsilon_) && ((y - x) < epsilon_);
}

bool DoubleEquality::operator()(const Complex& x, const Complex& y)
{
    return x == y;
}

