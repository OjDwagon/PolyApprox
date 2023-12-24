#include "complex.h"
#include "complexPolynomial.h"
#include "realPolynomial.h"
#include <chrono>

using namespace std;

#ifdef COMPLEX
    typedef ComplexPolynomial Polynomial;
    typedef Complex T;
#else
    typedef RealPolynomial Polynomial;
    typedef double T;
#endif


// Returns the average time over 50 runs of squaring the polynomial p
double multiplicationTimeTest(const Polynomial& p1, const Polynomial& p2, Polynomial (*multiplicationMethod)(const Polynomial&, const Polynomial&));

int main()
{
    vector<double> distributiveTimes, fourierTimes;
    vector<int> polynomialSizes({1, 2, 4, 8, 16, 32, 64, 128, 256, 194});

    for(int size: polynomialSizes){
        distributiveTimes.push_back(multiplicationTimeTest(Polynomial(vector<T>(size, 1)), Polynomial(vector<T>(size, 1)),  Polynomial::distributiveMultiplication));
        fourierTimes.push_back(multiplicationTimeTest(Polynomial(vector<T>(size, 1)), Polynomial(vector<T>(size, 1)), Polynomial::discreteFourierTransformMultiplication));
    }

    for(unsigned int i = 0; i < polynomialSizes.size(); ++i){
        cout << "Size: " << polynomialSizes[i] << " D: " << distributiveTimes[i] << " F: " << fourierTimes[i] << endl;
    }
    
    return 0;
}


double multiplicationTimeTest(const Polynomial& p1, const Polynomial& p2, Polynomial (*multiplicationMethod)(const Polynomial&, const Polynomial&))
{
    double total = 0;

    // Run 5 before beginning measurements to avoid cold start
    for(int i = 0; i < 5; ++i){
        multiplicationMethod(p1, p2);
    }

    for(int i = 0; i < 50; ++i){
        chrono::time_point start = chrono::steady_clock::now();
        multiplicationMethod(p1, p2);
        chrono::time_point end = chrono::steady_clock::now();

        chrono::duration<double, std::milli> duration = end - start;
        total += duration.count();
    }


    return total / 50;
}



