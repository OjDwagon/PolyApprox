
#include <vector>
#include <cassert>

#include "complex.h"
#include "complexPolynomial.h"
#include "complexPolynomial-test.h"

using namespace std;

typedef ComplexPolynomial Polynomial;

int main()
{
    evaluationTests();
    comparisonTests();
    additionTests();
    subtractionTests();
    basicMultiplicationTests();
    distributiveMultiplicationTests();
    fourierMultiplicationTests();
    derivativeTests();
    exponentiationTests();
    

    cout << "All tests passed." << endl;
    return 0;
}

void evaluationTests()
{
    // Real number evaluations
    Polynomial p1({1, 2, 3, 4, 5});

    assert(p1.evaluate(1) == 15);
    assert(p1(1.0) == p1.evaluate(1));
    assert(p1(0.5) == 3.5625);

    // Complex evaluations
    assert(p1(Complex(1, -1)) == Complex(-25, -16));

}

void additionTests()
{
    Polynomial p1;
    Polynomial p2({1, 2, 3, 4, 5, 6});

    for(int i = 0; i < 5; ++i){
        p2 += p1;
        assert(p2 == Polynomial({1, 2, 3, 4, 5, 6}));
    }

    assert((p2 + 1) == Polynomial({2, 2, 3, 4, 5, 6}));
    p2 += 1.0;
    assert(p2 == Polynomial({2, 2, 3, 4, 5, 6}));

    p2 += p2;
    assert(p2 == Polynomial({4, 4, 6, 8, 10, 12}));
}

void subtractionTests()
{
    Polynomial p1;
    Polynomial p2({1, 2, 3, 4, 5, 6});
    for(int i = 0; i < 5; ++i){
        p2 -= p1;
        assert(p2 == Polynomial({1, 2, 3, 4, 5, 6}));
    }

    assert((p2 - 1) == Polynomial({0, 2, 3, 4, 5, 6}));
    p2 -= 1.0;
    assert(p2 == Polynomial({0, 2, 3, 4, 5, 6}));

    p2 -= p2;
    assert(p2 == 0);
}

void basicMultiplicationTests()
{
    Polynomial p1({1, 2, 3, 4, 5});
    Polynomial p2(p1);

    assert((p1 * 1) == p1);
    assert((p1 * 2.0) == Polynomial({2, 4, 6, 8, 10}));
    
    p2 *= Complex(0, 1);
    assert(p2 == Polynomial({Complex(0, 1), Complex(0, 2), Complex(0, 3), Complex(0, 4), Complex(0, 5)}));
    p2 *= Complex(0, 1);
    assert(p2 == Polynomial({-1, -2, -3, -4, -5}));

    assert((0 * p1) == 0);
}

void distributiveMultiplicationTests()
{
    Polynomial p1({1, 2, 3, 4, 5});
    assert(Polynomial::distributiveMultiplication(p1, p1) == Polynomial({1, 4, 10, 20, 35, 44, 46, 40, 25}));

    Polynomial p2({Complex(0, 1), Complex(1, 3)});
    assert(Polynomial::distributiveMultiplication(p2, p2) == Polynomial({Complex(-1, 0), Complex(-6, 2), Complex(-8, 6)}));
}

void fourierMultiplicationTests()
{
    Polynomial p1({1, 2, 3, 4, 5});
    assert(Polynomial::discreteFourierTransformMultiplication(p1, p1) == Polynomial({1, 4, 10, 20, 35, 44, 46, 40, 25}));
    assert(Polynomial::discreteFourierTransformMultiplication(Polynomial({1, 1}), Polynomial({1, 1})) == Polynomial({1, 2, 1}));
    assert(Polynomial::discreteFourierTransformMultiplication(Polynomial({1, 3, 4}), Polynomial({Complex(), Complex(1, -1), Complex(3, 5)})) ==
            Polynomial({Complex(), Complex(1, -1), Complex(6, 2), Complex(13, 11), Complex(12, 20)}));
}

void derivativeTests()
{
    Polynomial p1({1, 2, 3, 4, 5});

    assert(p1.getDerivative() == Polynomial({2, 6, 12, 20}));
}

void comparisonTests()
{
    Polynomial p1({1, 2, 3, 4});
    Polynomial p2({1, 2, 3});

    assert(p1 != p2);

    p2.setCoefficient(3, 4);
    assert(p1 == p2);
}

void exponentiationTests()
{
    Polynomial p1({1, 1});

    assert(p1.exponentiate(5) == Polynomial({1, 5, 10, 10, 5, 1}));
    assert(p1.exponentiate(6) == Polynomial({1, 6, 15, 20, 15, 6, 1}));

    assert(Polynomial({-1, 1}).exponentiate(10) == Polynomial({1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1}));
}
