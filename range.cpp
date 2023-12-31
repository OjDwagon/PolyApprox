#include <iostream>
#include "range.h"

using namespace std;

ostream& operator<<(std::ostream& mycout, const Range& r)
{
    if(r.rangeType == Range::LEFT_INFINITE){
        cout << "A left-infinite range ";
    }
    else if(r.rangeType == Range::RIGHT_INFINITE){
        cout << "A right-infinite range ";
    }
    else{
        cout << "A finite range ";
    }

    if(r.increasing){
        cout << "increasing from: ";
    }
    else{
        cout << "decreasing from: ";
    }
    mycout << "(" << r.startX << ", " << r.startY << ") to (" << r.endX << ", " << r.endY << ")";
    
    return mycout;
}

Range::Range(double startIn, double startVal, double endIn, double endVal, double derivativeSample, int8_t typeOfRange):
    startX(startIn), startY(startVal), endX(endIn), endY(endVal), rangeType(typeOfRange)
{
    if(derivativeSample > 0){
        increasing = true;
    }
    else{
        increasing = false;
    }
}

Range::operator bool()
{   
    if(abs(startY) < 0.000001 || abs(endY) == 0.000001){
        return true;
    }
    else if(rangeType == LEFT_INFINITE){ // Left most range to infinity
        if(!increasing && endY < 0){ // If decreasing and ends below the x axis
            return true;
        }
        else if(endY > 0){ // If increasing and ends above the x axis
            return true;
        }
        else{
            return false;
        }
    }
    else if(rangeType == RIGHT_INFINITE){ // Right most range to infinity
        if(!increasing && startY > 0){ // If decreasing and starts above the x axis
            return true;
        }
        else if(startY < 0){ // If increasing and starts below the x axis
            return true;
        }
        else{
            return false;
        }
    }
    else{ // Finite range case
        if(!increasing && startY > 0 && endY < 0){
            return true;
        }
        else if(startY < 0 && endY > 0){
            return true;
        }
        else{
            return false;
        }
    }
   
}

const int8_t Range::LEFT_INFINITE = -1;
const int8_t Range::RIGHT_INFINITE = 1;
const int8_t Range::FINITE = 0;
