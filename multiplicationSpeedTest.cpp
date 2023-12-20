
#include "complex.h"
#include "complexPolynomial.h"
#include <chrono>

using namespace std;

// Returns the average time over 25 runs of squaring the polynomial p
double multiplicationTimeTest(const ComplexPolynomial& p, ComplexPolynomial (*multiplicationMethod)(const ComplexPolynomial&, const ComplexPolynomial&));

int main()
{
    vector<double> distributiveTimes, fourierTimes;
    vector<int> polynomialSizes({1, 2, 4, 8, 16, 32, 64, 128, 256, 192, 1000});

    for(int size: polynomialSizes){
        distributiveTimes.push_back(multiplicationTimeTest(ComplexPolynomial(vector<Complex>(size, 1)), ComplexPolynomial::distributiveMultiplication));
        fourierTimes.push_back(multiplicationTimeTest(ComplexPolynomial(vector<Complex>(size, 1)), ComplexPolynomial::discreteFourierTransformMultiplication));
    }

    for(unsigned int i = 0; i < polynomialSizes.size(); ++i){
        cout << "Size: " << polynomialSizes[i] << " D: " << distributiveTimes[i] << " F: " << fourierTimes[i] << endl;
    }
    
    return 0;
}


double multiplicationTimeTest(const ComplexPolynomial& p, ComplexPolynomial (*multiplicationMethod)(const ComplexPolynomial&, const ComplexPolynomial&))
{
    double total = 0;

    // Run 5 before beginning measurements to avoid cold start
    for(int i = 0; i < 5; ++i){
        multiplicationMethod(p, p);
    }

    for(int i = 0; i < 50; ++i){
        chrono::time_point start = chrono::steady_clock::now();
        multiplicationMethod(p, p);
        chrono::time_point end = chrono::steady_clock::now();

        chrono::duration<double, std::milli> duration = end - start;
        total += duration.count();
    }


    return total / 50;
}



