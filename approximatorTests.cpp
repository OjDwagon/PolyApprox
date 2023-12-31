
#include "approximator.h"
#include "realPolynomial.h"
#include <cassert>

using namespace std;

template<typename ChildType>
void approximatorTest(const GenericLinearApproximator<ChildType>& approximator);

int main()
{
    approximatorTest(IterativeLinearApproximator(150));
    approximatorTest(ErrorLinearApproximator(0.0001));

    cout << "All tests passed.";
    return 0;
}

template<typename ChildType>
void approximatorTest(const GenericLinearApproximator<ChildType>& approximator)
{
    // Quadratic estimation
    // cout << RealPolynomial({0, 0, 1}) << endl;
    assert(approximator(RealPolynomial({1, 0, 1})).size() == 0);
    assert(RealPolynomial(approximator(RealPolynomial({0, 0, 1}))) == RealPolynomial()); // y = x^2
    assert(RealPolynomial(approximator(RealPolynomial({-1, 0, 1}))) == RealPolynomial({-1, 1}));

    // Cubic Estimation
    assert(RealPolynomial(approximator(RealPolynomial({0, 0, 0, 1}))) == RealPolynomial()); // y = x^3
    assert(RealPolynomial(approximator(RealPolynomial({0, -5, 0, 1}))) == RealPolynomial({-2.236067977499, 0, 2.236067977499}));
    assert(RealPolynomial(approximator(RealPolynomial({6, 3, 2, 4}))) == RealPolynomial({-1.08424039378}));


    // 10th degree estimation
    assert(RealPolynomial(approximator(RealPolynomial({0, 0, 0, 0, -48, 20, 60, -23, -13, 3, 1}))) == RealPolynomial({-4, -3, -1, 0, 1, 2}));
}
