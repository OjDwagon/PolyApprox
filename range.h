#ifndef RANGE_H
#define RANGE_H

#include <iostream>

/**
 * @brief a Range struct that models monotonic intervals, does not model derivative intervals with derivative 0
*/
struct Range
{
    // enum finiteness {LEFT_INFINITE = -1, FINITE, RIGHT_INFINITE};
    static const int8_t LEFT_INFINITE;
    static const int8_t RIGHT_INFINITE;
    static const int8_t FINITE;

    double startX;
    double startY;
    double endX;
    double endY;
    bool increasing;
    int8_t rangeType; // -1 for infinite left range, 1 for infinite right range, 0 for finite range
    // int infinityType; // -1 for negative infinity, 1 for positive infinity, 0 for a finite number
    Range(double startIn, double startVal, double endIn, double endVal, double derivativeSample, int8_t typeOfRange = FINITE);
    operator bool(); // Returns true if a range crosses the x-axis (contains a root), assumes ranges are monotonic
};

/**
 * Prints out a Range to a given ostream
 * 
 * @relatesalso Range
*/
std::ostream& operator<<(std::ostream& mycout, const Range& r);

#endif
