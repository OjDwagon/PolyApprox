#define _USE_MATH_DEFINES
#include "complex.h"
#include "complex-test.h"
#include <cmath>
#include <cassert>

using namespace std;

bool doubleComp(double a, double b)
{
    return ((a - b) < 0.000001) && ((b - a) < 0.000001);
}

/**
 * Verifies a complex number by checking its charactarization against an expected value
*/
void verifyComplex(const Complex& num, const double& expectedX, const double& expectedY, const double& expectedR, const double& expectedTheta,
                    const bool& expectedCartesianStale, const bool& expectedPolarStale)
{
    // Get the actual values
    double actualX, actualY, actualR, actualTheta;
    bool actualCartesianStale, actualPolarStale;
    num.charactarization(actualX, actualY, actualR, actualTheta, actualCartesianStale, actualPolarStale);

    // Assert that they should be the same
    assert((actualX == expectedX) && (actualY == expectedY) && (actualCartesianStale == expectedCartesianStale));
    assert((actualX == expectedX) && (actualY == expectedY) && (actualPolarStale == expectedPolarStale));
}

int main(int argc, char* argv[])
{
    equalityTests();
    conversionTests();
    additionTests();
    subtractionTests();
    multiplicationTests();
    divisionTests();
    exponentiationTests();
    doublesCastingTests();

    cout << "All tests passed.";

    return 0;

}

void conversionTests()
{
    // Cartesian single coord constructor
    Complex num1(1, 1);
    Complex num2(sqrt(2), M_PI_4 + 2 * M_PI, false);

    assert(doubleComp(num1.getPolarR(), sqrt(2)));
    assert(doubleComp(num1.getPolarTheta(), M_PI_4));
    assert(num1 == num2);

    assert(Complex(0, 0) == Complex(0, 0, false));
}

void additionTests()
{
    Complex num1(1, 0);
    Complex num2(0, 1);
    Complex num3(1, 1);

    assert((num1 + num1) == Complex(2, 0));
    assert((num2 + num2) == Complex(0, 2));
    assert((num1 + num2) == num3);
    assert((num2 + 1) == num3);
    assert((1 + num2) == num3);


    for(int i = 2; i < 10; ++i){
        num3 += Complex(1, 1);
        assert(num3 == Complex(i, i));
    }

    num3 += 1;
    assert(num3 == Complex(10, 9));
}

void subtractionTests()
{
    Complex num1(1, 0);
    Complex num2(0, 1);
    Complex num3(1, 1);

    assert((num3 - num1) == num2);
    assert((num3 - num2) == num1);
    assert((1 - num3) == Complex(0, -1));

    for(int i = 0; i > -8; --i){
        num3 -= Complex(1, 1);
        assert(num3 == Complex(i, i));
    }

    num3 -= -2;
    assert(num3 == Complex(-5, -7));
}

void multiplicationTests()
{
    // Real number multiplication
    Complex num1(2, 0);
    Complex num2(-3, 0);
    assert(num1 * num1 == 4);
    assert(num2 * num2 == 9);
    assert(num1 * num2 == -6);
    assert(-5 * num1 == -10);

    // Complex multiplication
    Complex num3(0, 1);

    num3 *= Complex(0, 1);
    assert(num3 == -1);

    num3 *= Complex(0, 1);
    assert(num3 == Complex(0, -1));

    num3 *= Complex(0, 1);
    assert(num3 == 1);

    num3 *= Complex(0, 1);
    assert(num3 == Complex(0, 1));

    assert(2 * num3 == Complex(0, 2));

    assert(Complex(2, 3) * Complex(5, 4) == Complex(-2, 23));
}

void divisionTests()
{
    // Real number division
    Complex num1(2, 0);
    Complex num2(3, 0);

    assert((num1 / num2) == (double) 2 / 3);
    num2 /= -1;
    assert(num2 == -3);

    assert((num1 / num2) == ((double)-2 / 3));

    // Complex number division
    Complex num3(0, 1);

    assert((num3 / num3) == 1);
    assert((1 / num3) == Complex(0, -1));

    Complex num4(2, 5);
    Complex num5(6, 3);
    assert((num4 / num5) == (((double) 1 / 15) * Complex(9, 8)));
}

void equalityTests()
{
    // Complex with complex
    Complex num1(1.0, 3.1);
    Complex num2(1.0, 3.1);
    assert(num1 == num2); // 1 + 3.1i == 2 + 3.1i

    Complex num3(0, 5);
    Complex num4(5, M_PI_2, false);
    assert(num3 == num4); // Same number, one from cartesian, one from polar

    assert(num1 != num3); // Different numbers

    // Complex with double
    Complex num5(1.0, 1.0);
    Complex num6(1.0, 0);
    
    assert(num5 != 1); // 1 + i != 1
    assert(num6 == 1); // 1 == 1
}

void exponentiationTests()
{
    // Real number exponentiation
    assert((Complex(2, 0) ^ 2.5) == pow(2, 2.5));
    assert((Complex(13, 5.5) ^ 0) == 1);

    // Complex exponentiation
    assert((Complex(1, -1) ^ 5) == Complex(-4, 4));

    Complex num1(0, 1);

    for(int i = 0; i < 10; ++i){ // i ^ 4 = 1
        assert((num1 ^ (4 * i)) == 1);
    }

    for(int i = 0; i < 10; ++i){ // i ^ 1 = i
        assert((num1 ^ (4 * i + 1)) == Complex(0, 1));
    }

    for(int i = 0; i < 10; ++i){ // i ^ 2 = -1
        assert((num1 ^ (4 * i + 2)) == -1);
    }

    for(int i = 0; i < 10; ++i){ // i ^ 3 = -i
        assert((num1 ^ (4 * i + 3)) == Complex(0, -1));
    }
}

void doublesCastingTests()
{
    Complex num1((double) 1 / 3, 1);
    assert(doubleComp(static_cast<double>(num1), (double) 1 / 3));
    assert(doubleComp(num1, (double) 1 / 3));
}


