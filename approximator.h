#ifndef APPROXIMATOR_H
#define APPROXIMATOR_H
#include <vector>
#include "realPolynomial.h"

/**
 * @brief An abstract linear root approximator that implements the curiously recurrent template pattern, expects child to implement an approximate method
 * @tparam ChildType the childtype with an approximate method
*/
template<typename ChildType>
class GenericLinearApproximator
{
    public:
        /**
         * Approximates all roots of a Polynomial using Linear Approximation by generating appropriate guesses
         * 
         * @param f the Polynomial whose roots will be approximated
         * @returns a vector of the approximated roots
        */
        std::vector<double> operator()(const RealPolynomial& f) const;

        /**
         * Approximates a single root of a Polynomial using Linear Approximation, beginning at a guess
         * 
         * @param f the Polynomial that will have its root approximated
         * @param guess the guess where linear approximation will begin from
         * @returns the approximated root
        */
        double operator()(const RealPolynomial& f, double& guess) const;

    protected:
        struct DoubleLess
        {
            bool operator()(double a, double b) const
            {
                return b - a > 0.000001;
            }
        };

        struct DoubleEqual
        {
            bool operator()(double a, double b) const
            {
                return ((a - b) < 0.000001) && ((b - a) < 0.000001);
            }
        };
};


template<typename ChildType>
std::vector<double> GenericLinearApproximator<ChildType>::operator()(const RealPolynomial& f) const
{
    return static_cast<const ChildType*>(this)->approximate(f);
}

template<typename ChildType>
double GenericLinearApproximator<ChildType>::operator()(const RealPolynomial& f, double& guess) const
{
    return static_cast<const ChildType*>(this)->approximate(f, guess);
}

/**
 * @brief An implementation of the generic linear root approximator which uses a fixed number of iterations to generate approximations
*/
class IterativeLinearApproximator: public GenericLinearApproximator<IterativeLinearApproximator>
{
    public:
        IterativeLinearApproximator(unsigned int i = 10);
        void setIterations(unsigned int i);
        unsigned int getIterations() const;
        std::vector<double> approximate(const RealPolynomial& f) const;
        double approximate(const RealPolynomial& f, double& guess) const;
    private:
        double approximate(const RealPolynomial& f, const RealPolynomial& fprime, double& guess, unsigned int iterations) const;
        unsigned int iterations_;
};

/**
 * @brief An implementation of the generic root linear approximator which continues using linear approximation until
 * approximations are within a specified error range from 0
*/
class ErrorLinearApproximator: public GenericLinearApproximator<ErrorLinearApproximator>
{
    public:
        ErrorLinearApproximator(double error = 0.001, unsigned int maxIterations = 100);
        void setError(double error);
        double getError();
        void setMaxIterations(unsigned int maxIterations);
        unsigned int getMaxIterations();
        std::vector<double> approximate(const RealPolynomial& f) const;
        double approximate(const RealPolynomial& f, double& guess) const;
    private:
        double approximate(const RealPolynomial& f, const RealPolynomial& fprime, double& guess, unsigned int iterations) const;
        double error_;
        unsigned int maxIterations_;

};

std::vector<double> getRoots(RealPolynomial f); // Returns the roots of a polynomial less than degree 2 in sorted order
double avg(double x, double y); // Returns the average of two numbers

#endif
