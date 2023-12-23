#include "ComplexPolynomial.h"
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

ComplexPolynomial::ComplexPolynomial():
    GenericPolynomial<Complex, ComplexPolynomial>()
{

}

ComplexPolynomial::ComplexPolynomial(const vector<Complex>& coefficients):
    GenericPolynomial<Complex, ComplexPolynomial>(coefficients)
{

}

ComplexPolynomial::ComplexPolynomial(const initializer_list<Complex>& coefficients):
    GenericPolynomial<Complex, ComplexPolynomial>(coefficients)
{

}

ComplexPolynomial ComplexPolynomial::operator*(const ComplexPolynomial& rhs) const
{
    if(max(degree_, rhs.degree_) > 190){ // Discrete Fourier transform is more efficient once multiplying over 190 degrees
        return discreteFourierTransformMultiplication(*this, rhs);
    }
    else{
        return distributiveMultiplication(*this, rhs);
    }
}

ComplexPolynomial ComplexPolynomial::operator*(Complex rhs) const
{
    return *(static_cast<const GenericPolynomial<Complex, ComplexPolynomial>*>(this)) * rhs;
}

ostream& operator<<(ostream& out, const ComplexPolynomial& rhs)
{
    out << ComplexPolynomial::getPrintableComplex(rhs.coefficients_[0]);
    for(unsigned int i = 1; i <= rhs.degree_; ++i){
        out << " + " << ComplexPolynomial::getPrintableComplex(rhs.coefficients_[i]) << "x^" << i; 
    }
    return out;
}

string ComplexPolynomial::getPrintableComplex(const Complex& num)
{
    stringstream ss;
    ss << "(" << num.getCartesianX() << " + " << num.getCartesianY() << "i)";
    return ss.str();
}

ComplexPolynomial ComplexPolynomial::discreteFourierTransformMultiplication(const ComplexPolynomial& inLhs, const ComplexPolynomial& inRhs)
{   
    ComplexPolynomial lhs(inLhs);
    ComplexPolynomial rhs(inRhs);

    // Assumes lhs is a higher degree, pad up to the next power of 2 for even recursive calls
    unsigned int n = getNextPowerOfTwo(lhs.degree_ + 1); // n is the number of coefficients in each function
    lhs.pad(n - 1);
    rhs.pad(n - 1);

    // Precompute the needed roots of unity
    vector<vector<Complex>> rootsOfUnity(log2(n) + 1); // A vector where the ith entry is the 2ith roots of unity
    for(unsigned int i = 1, j = 0; i <= n; i *= 2, ++j){
        rootsOfUnity[j] = calcRootsOfUnity(2 * i);
    }

    // Recursively evaluate the lhs and rhs on the 2nth roots of unity
    vector<Complex> lhsAtRoots = evaluateAtRootsOfUnity(lhs, 2 * n, rootsOfUnity);
    vector<Complex> rhsAtRoots = evaluateAtRootsOfUnity(rhs, 2 * n, rootsOfUnity);

    // Use this to evaluate the product on the roots of unity
    vector<Complex> productAtRoots(2 * n);
    for(unsigned int i = 0; i < 2 * n; ++i){
        productAtRoots[i] = lhsAtRoots[i] * rhsAtRoots[i];
    }

    // Use an intermediate fucntion to generate the final coefficients
    ComplexPolynomial intermediateFunction(productAtRoots);
    intermediateFunction.pad(2 * n - 1);
    vector<Complex> intermediateRoots = evaluateAtRootsOfUnity(intermediateFunction, 2 * n, rootsOfUnity);
    vector<Complex> finalCoefficients(2 * n);

    finalCoefficients[0] = intermediateRoots[0] / ((double)finalCoefficients.size());

    for(unsigned int i = 1; i < 2 * n; ++i){
        finalCoefficients[i] = intermediateRoots[2 * n - i] / ((double)finalCoefficients.size());
    }

    return ComplexPolynomial(finalCoefficients);
}

vector<Complex> ComplexPolynomial::evaluateAtRootsOfUnity(const ComplexPolynomial& f, unsigned int n, const vector<vector<Complex>>& rootsOfUnity)
{
    vector<Complex> answer(n);
    const vector<Complex>& currentRoots = rootsOfUnity[log2(n) - 1];

    // Base Case
    if(n < 3){
        for(unsigned int i = 0; i < n; ++i){
            answer[i] = f(currentRoots[i]);
        }
        return answer;
    }

    // Recursive Step

    // Generate even and odd coefficient functions
    vector<Complex> evenCoefficients;
    vector<Complex> oddCoefficients;

    // Guaranteed to have these 2 coefficients for initialization
    evenCoefficients.push_back(f.getCoefficient(0));
    oddCoefficients.push_back(f.getCoefficient(1));

    // Assign the rest of the coefficients
    for(unsigned int i = 2; i <= f.getDegree(); ++i){
        if(i % 2 == 0){
            evenCoefficients.push_back(f.getCoefficient(i));
        }
        else{
            oddCoefficients.push_back(f.getCoefficient(i));
        }
    }

    ComplexPolynomial even(evenCoefficients);
    ComplexPolynomial odd(oddCoefficients);
    even.pad((f.getDegree() - 1) / 2);
    odd.pad((f.getDegree() - 1) / 2);

    // Recursively evaluate evens and odds on nth roots of unity
    vector<Complex> evenAtRoots = evaluateAtRootsOfUnity(even, n / 2, rootsOfUnity);
    vector<Complex> oddAtRoots = evaluateAtRootsOfUnity(odd, n / 2, rootsOfUnity);

    // Combine answers to form final
    for(unsigned int i = 0; i < n; ++i){
        answer[i] = evenAtRoots[i % (n / 2)] + currentRoots[i] * oddAtRoots[i % (n / 2)];
    }
    
    return answer;
}


vector<Complex> ComplexPolynomial::calcRootsOfUnity(unsigned int n)
{
    vector<Complex> roots;

    for(unsigned int i = 0; i < n; ++i){
        roots.push_back(Complex(1, (2 * M_PI * i) / n, false));
    }

    return roots;
}

unsigned int ComplexPolynomial::getNextPowerOfTwo(unsigned int n)
{
    return static_cast<unsigned int>(pow(2, ceil(log2(n))));
}


ComplexPolynomial operator+(Complex lhs, const ComplexPolynomial& rhs)
{
    return rhs + lhs;
}

ComplexPolynomial operator-(Complex lhs, const ComplexPolynomial& rhs)
{
    return -1 * rhs + lhs;
}

ComplexPolynomial operator*(Complex lhs, const ComplexPolynomial& rhs)
{
    return static_cast<GenericPolynomial<Complex, ComplexPolynomial>>(rhs) * lhs;
}


bool operator==(Complex lhs, const ComplexPolynomial& rhs)
{
    return rhs == lhs;
}

bool operator!=(Complex lhs, const ComplexPolynomial& rhs)
{
    return rhs != lhs;
}
