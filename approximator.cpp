#include "approximator.h"
#include "realPolynomial.h"
#include "range.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <set>

using namespace std;

IterativeLinearApproximator::IterativeLinearApproximator(unsigned int i):
    iterations_(i)
{

}

void IterativeLinearApproximator::setIterations(unsigned int i)
{
    iterations_ = i;
}

unsigned int IterativeLinearApproximator::getIterations() const
{
    return iterations_;
}

vector<double> IterativeLinearApproximator::approximate(const RealPolynomial& f) const
{
    if(f.getDegree() < 3){ // If a closed form exists return roots immediately;
        return getRoots(f);
    }
    
    // Generate up to degree - 1 guesses
    vector<double> guesses;
    RealPolynomial fprime = f.getDerivative(); // Get the derivative
    vector<double> critValues = approximate(fprime); // Find crit values using the same number of iterations

    // Generate ranges based on the critical values
    double a = abs(f.getCoefficient(f.getDegree())); // Scaling factor based on coefficient of highest order term (used for infinity ranges)
    vector<Range> ranges;

    if(critValues.size() == 0){ // No critical values means it must be a line
        ranges.push_back(Range(-1 * a, f(-1 * a), a, f(a), fprime(0))); // Use a offsets for infinity
    }
    else{
        // Set the beginning range
        double leftMostCrit = critValues[0];
        double rightMostCrit = critValues[critValues.size() - 1];
        ranges.push_back(Range(leftMostCrit - a, f(leftMostCrit - a), leftMostCrit, f(leftMostCrit), fprime(avg(leftMostCrit, leftMostCrit - a)), Range::LEFT_INFINITE));

        // Set the rest of the ranges except the last
        for(unsigned int i = 1; i < critValues.size(); ++i){
            ranges.push_back(Range(critValues[i - 1], f(critValues[i - 1]), critValues[i], f(critValues[i]), fprime(avg(critValues[i - 1], critValues[i]))));
        }

        // Set the last range
        ranges.push_back(Range(rightMostCrit, f(rightMostCrit), 
                rightMostCrit + a, f(rightMostCrit + a), 
                fprime(avg(rightMostCrit, rightMostCrit + a)), Range::RIGHT_INFINITE));

    }

    // Iterate through all the ranges, check if they contain zero using the bool() operator overload
    // If the range happens to end at zero, record that as a root
    set<double, GenericLinearApproximator<IterativeLinearApproximator>::DoubleLess> rootSet;
    for(unsigned int i = 0; i < ranges.size(); ++i){
        if(ranges[i]){
            DoubleEqual comp;
            if(comp(ranges[i].startY, 0)){
                rootSet.insert(ranges[i].startX);
            }
            else if(comp(ranges[i].endY, 0)){
                rootSet.insert(ranges[i].endX);
            }
            else{
                guesses.push_back(avg(ranges[i].startX, ranges[i].endX)); // Add the center of the range as a guess
            }
        }
    }

    // Iterate through all guesses to produce root estimates
    for(unsigned int i = 0; i < guesses.size(); ++i){
        rootSet.insert(approximate(f, fprime, guesses[i], iterations_));
    }

    return vector<double>(rootSet.begin(), rootSet.end());
}

double IterativeLinearApproximator::approximate(const RealPolynomial& f, double& guess) const
{
    return approximate(f, f.getDerivative(), guess, iterations_);
}

double IterativeLinearApproximator::approximate(const RealPolynomial& f, const RealPolynomial& fprime, double& guess, unsigned int iterations) const
{
    if(iterations == 0){ // Base case
        return guess;
    }
    if(fprime(guess) == 0){
        throw domain_error("Can't perform linear approximation when derivative is 0");
    }

    double nextGuess = -1 * (f(guess) / fprime(guess)) + guess; // Calculate the next step and recurse
    return approximate(f, fprime, nextGuess, iterations - 1);
}

ErrorLinearApproximator::ErrorLinearApproximator(double error, unsigned int maxIterations):
    error_(error), maxIterations_(maxIterations)
{
    if(error_ < 0){ // Error must be positive
        error_ *= -1;
    }
}

void ErrorLinearApproximator::setError(double error)
{
    error_ = error;
    if(error_ < 0){ // Error must be positive
        error_ *= -1;
    }
}

double ErrorLinearApproximator::getError()
{
    return error_;
}

void ErrorLinearApproximator::setMaxIterations(unsigned int maxIterations)
{
    maxIterations_ = maxIterations;
}

unsigned int ErrorLinearApproximator::getMaxIterations()
{
    return maxIterations_;
}

vector<double> ErrorLinearApproximator::approximate(const RealPolynomial& f) const
{
    if(f.getDegree() < 3){ // If a closed form exists return roots immediately;
        return getRoots(f);
    }
    
    // Generate up to degree -1 guesses
    vector<double> guesses;
    RealPolynomial fprime = f.getDerivative(); // Get the derivative
    vector<double> critValues = approximate(fprime); // Find crit values using the same number of iterations

    // Generate ranges based on the critical values
    double a = abs(f.getCoefficient(f.getDegree())); // Scaling factor based on coefficient of highest order term (used for infinity ranges)
    vector<Range> ranges;

    if(critValues.size() == 0){ // No critical values means it must be a line
        ranges.push_back(Range(-1 * a, f(-1 * a), a, f(a), fprime(0))); // Default to 10 offsets for infinity
    }
    else{
        // Set the beginning range
        double leftMostCrit = critValues[0];
        double rightMostCrit = critValues[critValues.size() - 1];
        ranges.push_back(Range(leftMostCrit - a, f(leftMostCrit - a), leftMostCrit, f(leftMostCrit), fprime(avg(leftMostCrit, leftMostCrit - a)), Range::LEFT_INFINITE));

        // Set the rest of the ranges except the last
        for(unsigned int i = 1; i < critValues.size(); ++i){
            ranges.push_back(Range(critValues[i - 1], f(critValues[i - 1]), critValues[i], f(critValues[i]), fprime(avg(critValues[i - 1], critValues[i]))));
        }

        // Set the last range
        ranges.push_back(Range(rightMostCrit, f(rightMostCrit), 
                rightMostCrit + a, f(rightMostCrit + a), 
                fprime(avg(rightMostCrit, rightMostCrit + a)), Range::RIGHT_INFINITE));

    }


    // Iterate through all the ranges, check if they contain zero using the bool() operator overload
    // If the range happens to end at zero, record that as a root
    set<double, GenericLinearApproximator<ErrorLinearApproximator>::DoubleLess> rootSet;
    for(unsigned int i = 0; i < ranges.size(); ++i){
        if(ranges[i]){
            DoubleEqual comp;
            if(comp(ranges[i].startY, 0)){
                rootSet.insert(ranges[i].startX);
            }
            else if(comp(ranges[i].endY, 0)){
                rootSet.insert(ranges[i].endX);
            }
            else{
                guesses.push_back(avg(ranges[i].startX, ranges[i].endX)); // Add the center of the range as a guess
            }
        }
    }

    // Iterate through all guesses to produce root estimates
    for(unsigned int i = 0; i < guesses.size(); ++i){
        rootSet.insert(approximate(f, fprime, guesses[i], maxIterations_));
    }

    return vector<double>(rootSet.begin(), rootSet.end());
}

double ErrorLinearApproximator::approximate(const RealPolynomial& f, double& guess) const
{
    return approximate(f, f.getDerivative(), guess, maxIterations_);
}

double ErrorLinearApproximator::approximate(const RealPolynomial& f, const RealPolynomial& fprime, double& guess, unsigned int iterations) const
{
    if(iterations == 0 || (abs(f(guess)) < error_)){ // Base case (check if error is low enough or if we have gone through the max number of iterations);
        return guess;
    }
    if(fprime(guess) == 0){
        throw domain_error("Can't perform linear approximation when derivative is 0");
    }

    double nextGuess = -1 * (f(guess) / fprime(guess)) + guess; // Calculate the next step and recurse
    return approximate(f, fprime, nextGuess, iterations - 1);
}

vector<double> getRoots(RealPolynomial f)
{
    vector<double> roots;
    switch(f.getDegree()){
        case 0:
            return roots;
        case 1: // Linear
            roots.push_back((-1 * f.getCoefficient(0)) / f.getCoefficient(1));
            return roots;
        case 2: // Quadratic
            double a = f.getCoefficient(2);
            double b = f.getCoefficient(1);
            double c = f.getCoefficient(0);
            double discriminant = b * b - 4 * a * c;
            if(discriminant > 0){ // Check discriminant
                roots.push_back((-1 * b + sqrt(b * b - 4 * a * c)) / (2 * a));
                roots.push_back((-1 * b - sqrt(b * b - 4 * a * c)) / (2 * a));
                if(abs(roots[0] - roots[1]) < 0.000001){ // Check if the roots are duplicated, uses double comparison
                    roots.pop_back();
                }
            }
            else if(discriminant == 0){
                roots.push_back(b / (2 * a));
            }
            sort(roots.begin(), roots.end());
            return roots;
    }
    // I'm considering 3rd degree polynomial too complicated for now
    throw logic_error("No easy closed form to get the roots of the function.");
}

double avg(double x, double y){
    return (x + y) / 2;
}
