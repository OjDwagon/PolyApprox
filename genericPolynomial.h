#ifndef GENERIC_POLYNOMIAL_H
#define GENERIC_POLYNOMIAL_H

#include <vector>

/**
 * @brief A generic polynomial class modeling polynomials of non-negative degree and their operators (excluding division)
 * Designed to be inherited from to support different types such as complex numbers
 * Children must implement an TrueType operator*(const TrueType&) overload and their own constructors
 * Children also must implement any required friend functions and global operator overloads
 * 
 * @tparam T, the type of the coefficients of the polynomial, must have the following defined:
 *          default constructor T() which evaluates to the 0 element 
 *          any number T^0 = T(1), the multiplicative identity element
 *          addition, subtraction, multiplication
 * 
 * @tparam TrueType, the type of the inheriting child, used to ensure methods return the appropriate child type
 *          children must implement the TrueType operator*(const TrueType&) method
 *          to implement TrueType operator*(const T&) method, children can cast this to be a parent and use the parent implementation
*/
template <typename T, typename TrueType>
class GenericPolynomial
{
    public:
        /**
         * Default constructor for a TrueType, defaults to the constant 0
        */
        GenericPolynomial();

        /**
         * GenericPolynomial constructor given coefficients
         * 
         * @tparam V the type of the coefficients, must support casting to a Complex number
         * @param coefficients a vector of coefficeints, the ith element of the vector is used as the coefficient for x^i
        */
        template <typename V>
        GenericPolynomial(const std::vector<V>& coefficients);

        /**
         * GenericPolynomial constructor given coefficients
         * 
         * @tparam V the type of the coefficients, must support casting to a Complex number
         * @param coefficients an initinialization of complex numbers, the ith element of the vector is used as the coefficient for x^i
        */
        template <typename V>
        GenericPolynomial(const std::initializer_list<V>& coefficients);       

        /**
         * Sets the coefficient corresponding to x^power to be val, does nothing if a negative power is given
         * 
         * @param power the degree of x corresponding to the coefficient to be set
         * @param val the value that the coefficient will be set to
        */
        void setCoefficient(unsigned int power, T val); // Does nothing if an invalid power is given

        /**
         * A linear time evaluation method of the polynomial
         * 
         * @tparam V the type of the input into the polynomial, must have defined arithmetic operations with Complex Numbers
         * @param x the value that the polynomial will be evaluated at
         * 
         * @returns the value of the polynomial evaluated at x as a complex number
        */
        template <typename V>
        T evaluate(V x) const;

        /**
         * Operator() that allows for evaluation of a polynomial in function style syntax
         * 
         * @tparam V the type of the input into the polynomial, must have defined arithmetic operations with Complex Numbers
         * @param x the value that the polynomial will be evaluated at
         * 
         * @returns the value of the polynomial evaluated at x as a complex number
        */
        template <typename V>
        T operator()(V x) const;

        TrueType operator+(const TrueType& rhs) const; // Adds 2 polynomials
        TrueType operator+(const T& rhs) const;
        TrueType& operator+=(const TrueType& rhs);
        TrueType& operator+=(const T& rhs);    

        TrueType operator-(const TrueType& rhs) const;
        TrueType operator-(const T& rhs) const;        
        TrueType& operator-=(const TrueType& rhs);
        TrueType& operator-=(const T& rhs);     

        TrueType operator*(const TrueType& rhs) const; // This operator must be implemented in children
        TrueType operator*(const T& rhs) const; // This operator must be implemented in children 
        TrueType& operator*=(const TrueType& rhs);
        TrueType& operator*=(const T& rhs);

        bool operator==(const TrueType& rhs) const;
        bool operator!=(const TrueType& rhs) const;
        bool operator==(double rhs) const;
        bool operator!=(double rhs) const;    

        /**
         * Calculates and returns the square of this polynomial
         * 
         * @returns a GenericPolynomial which is the result of squaring this polynomial
        */
        TrueType squared() const;

        /**
         * Calculates and returns the polynommial raised to the desired non-negative power
         * 
         * @param power a non-negative number which is the power the polynomial will be raised to
         * 
         * @returns a GenericPolynomial which is the result of raising this polynomial to the specified power
        */
        TrueType exponentiate(unsigned int power) const;

        /**
         * Getter for the coefficient associated with x^power, returns 0 if given an invalid coefficient
         * 
         * @param power the degree of x corresponding to the coefficient to be set
         * 
         * @returns the coefficient associated with x^power
        */
        T getCoefficient(unsigned int power) const;

        /**
         * Getter for the degree of the polynomial, which is the highest n such that the coefficient for x^n is not 0
         * 
         * @returns the degree of the GenericPolynomial
        */
        unsigned int getDegree() const;

        /**
         * Applies the power rule to calculate the derivative of this polynomial
         * 
         * @returns a GenericPolynomial which is the derivate of this polynomial
        */
        TrueType getDerivative() const; 

        /**
         * Returns the product of 2 polynomials, implemented using standard distribution in quadratic time
         * 
         * @param lhs the first polynomial to multiply
         * @param rhs the second polynomial to multiply
         * 
         * @returns the polynomial lhs * rhs
        */
        static TrueType distributiveMultiplication(const TrueType& lhs, const TrueType& rhs);

    protected:
        unsigned int degree_;
        std::vector<T> coefficients_;
        static double epsilon_; // Constant used for double comparison, doubles different by less than this amount are considered equal

        /**
         * Removes all leading zero coefficients (ignoring the constant coefficent x^0)
        */
        void unpad();

        /**
         * Adds leading zero coefficients to a polynomial to make it seem like a specified degree
         * 
         * @param degree the desired degree of the polynomial after the padding operation
        */
        void pad(unsigned int degree);

        /**
         * Helper function for double comparison, returns true if the doubels are within epsilon_ of each other
         * 
         * @param a, b the numbers to be compared
         * 
         * @returns true if a and b are within epsilon_ of each other
        */
        bool static doubleComp(double a, double b);
};

template<typename T, typename TrueType> template <typename V>
GenericPolynomial<T, TrueType>::GenericPolynomial(const std::vector<V>& coefficients):
    degree_(coefficients.size() - 1), coefficients_(coefficients.begin(), coefficients.end())
{
    unpad();
}

template<typename T, typename TrueType> template <typename V>
GenericPolynomial<T, TrueType>::GenericPolynomial(const std::initializer_list<V>& coefficients):
    degree_(coefficients.size() - 1)
{
    for(const V& c: coefficients){
        coefficients_.push_back(T(c));
    }

    unpad();
}

template<typename T, typename TrueType> template <typename V>
T GenericPolynomial<T, TrueType>::evaluate(V x) const
{
    T answer = T();
    for(int i = degree_; i > 0; --i){ // View the polynomial in a factored form
        answer += coefficients_[i];
        answer *= x;
    }

    answer += coefficients_[0];

    return answer;
}

template<typename T, typename TrueType> template <typename V>
T GenericPolynomial<T, TrueType>::operator()(V x) const
{
    return evaluate(x);
}

template<typename T, typename TrueType>
GenericPolynomial<T, TrueType>::GenericPolynomial():
    degree_(0), coefficients_({T()})
{

}

template<typename T, typename TrueType>
void GenericPolynomial<T, TrueType>::setCoefficient(unsigned int power, T val)
{
    if(power <= degree_){ // If it is a valid coefficient
        coefficients_[power] = val;
    }
    else{ // If it needs to be expanded 
        while(coefficients_.size() < power){
            coefficients_.push_back(T());
        }
        coefficients_.push_back(val);
        degree_ = power;
    }

    unpad();
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::operator+(const TrueType& rhs) const
{
    // Assumes lhs is bigger, if rhs is bigger flip them
    if(rhs.degree_ > degree_){
        return rhs + *(static_cast<const TrueType*>(this));
    }

    TrueType lhs = *(static_cast<const TrueType*>(this)); // Make a copy

    for(unsigned int i = 0; i <= rhs.degree_; ++i){ // Add all the coefficients
        lhs.coefficients_[i] += rhs.coefficients_[i];
    }

    return lhs; // Return it
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::operator+(const T& rhs) const
{
    TrueType lhs = *(static_cast<const TrueType*>(this)); // Make a copy
    lhs.coefficients_[0] += rhs; // Add the constant
    return lhs; // Return it
}

template<typename T, typename TrueType>
TrueType& GenericPolynomial<T, TrueType>::operator+=(const TrueType& rhs)
{
    if(rhs.degree_ > degree_){ // Pad up if necessary
        pad(rhs.degree_);
    }

    for(unsigned int i = 0; i <= rhs.degree_; ++i){ // Add all the coefficients
        coefficients_[i] += rhs.coefficients_[i];
    }

    unpad();

    return *(static_cast<TrueType*>(this));
}

template<typename T, typename TrueType>
TrueType& GenericPolynomial<T, TrueType>::operator+=(const T& rhs)
{
    coefficients_[0] += rhs; // Add the constant
    return *(static_cast<TrueType*>(this));
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::operator-(const TrueType& rhs) const
{
    TrueType lhs = *(static_cast<const TrueType*>(this)); // Make a copy

    if(rhs.degree_ > degree_){ // Pad up if necessary
        lhs.pad(rhs.degree_);
    }

    for(unsigned int i = 0; i <= rhs.degree_; ++i){ // Add all the coefficients
        lhs.coefficients_[i] -= rhs.coefficients_[i];
    }

    lhs.unpad();

    return lhs; // Return it
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::operator-(const T& rhs) const
{
    TrueType lhs = *(static_cast<const TrueType*>(this)); // Make a copy
    lhs.coefficients_[0] -= rhs; // Add the constant
    return lhs; // Return it
}

template<typename T, typename TrueType>
TrueType& GenericPolynomial<T, TrueType>::operator-=(const TrueType& rhs)
{
    if(rhs.degree_ > degree_){ // Pad up if necessary
        pad(rhs.degree_);
    }

    for(unsigned int i = 0; i <= rhs.degree_; ++i){ // Add all the coefficients
        coefficients_[i] -= rhs.coefficients_[i];
    }

    unpad();

    return *(static_cast<TrueType*>(this));
}

template<typename T, typename TrueType>
TrueType& GenericPolynomial<T, TrueType>::operator-=(const T& rhs)
{
    coefficients_[0] -= rhs; // Subtract the constant
    return *(static_cast<TrueType*>(this));
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::operator*(const TrueType& rhs) const
{
    return (*(static_cast<const TrueType*>(this))) * rhs;
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::operator*(const T& rhs) const
{
    std::vector<T> product = (static_cast<const TrueType*>(this))->coefficients_; // Make a copy

    for(unsigned int i = 0; i <= (static_cast<const TrueType*>(this))->degree_; ++i){ // Multiply everything
        product[i] *= rhs;
    }


    return TrueType(product); // Return
}

template<typename T, typename TrueType>
TrueType& GenericPolynomial<T, TrueType>::operator*=(const TrueType& rhs)
{
    TrueType product = *(static_cast<const TrueType*>(this)) * rhs;
    coefficients_ = product.coefficients_;
    degree_ = product.degree_;

    return *(static_cast<TrueType*>(this));
}

template<typename T, typename TrueType>
TrueType& GenericPolynomial<T, TrueType>::operator*=(const T& rhs)
{
    for(unsigned int i = 0; i <= degree_; ++i){
        coefficients_[i] *= rhs;
    }

    unpad();

    return *(static_cast<TrueType*>(this));
}


template<typename T, typename TrueType>
bool GenericPolynomial<T, TrueType>::operator==(const TrueType& rhs) const
{
    // Two polynomials of different degree can't be equal
    if(degree_ != rhs.degree_){
        return false;
    }

    // Check if all coefficients match
    for(unsigned int i = 0; i < coefficients_.size(); ++i){
        if(coefficients_[i] != rhs.coefficients_[i]){
            return false;
        }
    }

    return true;
    
}

template<typename T, typename TrueType>
bool GenericPolynomial<T, TrueType>::operator!=(const TrueType& rhs) const
{
    // Two polynomials of different degree can't be equal
    if(degree_ != rhs.degree_){
        return true;
    }

    // Check if all coefficients match
    for(unsigned int i = 0; i < coefficients_.size(); ++i){
        if(coefficients_[i] == rhs.coefficients_[i]){
            return false;
        }
    }

    return true;
}

template<typename T, typename TrueType>
bool GenericPolynomial<T, TrueType>::operator==(double rhs) const
{
    if(degree_ == 0){
        return rhs == coefficients_[0];
    }

    return false;
}

template<typename T, typename TrueType>
bool GenericPolynomial<T, TrueType>::operator!=(double rhs) const
{
    if(degree_ != 0){
        return true;
    }

    return rhs != coefficients_[0];
}

template<typename T, typename TrueType>
bool operator==(double lhs, const TrueType& rhs)
{
    return rhs == lhs;
}

template<typename T, typename TrueType>
bool operator!=(double lhs, const TrueType& rhs)
{
    return rhs != lhs;
}

template<typename T, typename TrueType>
bool operator==(int lhs, const TrueType& rhs)
{
    return rhs == lhs;
}

template<typename T, typename TrueType>
bool operator!=(int lhs, const TrueType& rhs)
{
    return rhs != lhs;
}



template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::squared() const
{
    return (*(static_cast<const TrueType*>(this))) * (*(static_cast<const TrueType*>(this)));
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::exponentiate(unsigned int power) const
{   
    // Recursively multiply
    if(power == 0){ // Base Case
        return TrueType(std::vector<T>({T(1)}));
    }
    else if(power == 1){ // Base Case
        return *(static_cast<const TrueType*>(this));
    }
    else if(power == 2){ // Base Case
        return (static_cast<const TrueType*>(this))->squared();
    }
    else if(power % 2 == 0){ // Even Recursive Case
        return (exponentiate(power / 2)).squared();
    }
    else{ // Odd Recursive case
        return *(static_cast<const TrueType*>(this)) * ((exponentiate(power / 2)).squared());
    }
}

template<typename T, typename TrueType>
T GenericPolynomial<T, TrueType>::getCoefficient(unsigned int power) const
{
    if(power <= degree_){
        return coefficients_[power];
    }
    else{
        return T();
    }
}

template<typename T, typename TrueType>
unsigned int GenericPolynomial<T, TrueType>::getDegree() const
{
    return degree_;
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::getDerivative() const
{
    std::vector<T> derivative(degree_);

    for(unsigned int i = 0; i < degree_; ++i){
        derivative[i] = (int)(i + 1) * coefficients_[i + 1];
    }

    return TrueType(derivative);
}

template<typename T, typename TrueType>
TrueType GenericPolynomial<T, TrueType>::distributiveMultiplication(const TrueType& lhs, const TrueType& rhs)
{
    TrueType product;
    product.pad(lhs.degree_ + rhs.degree_);

    for(unsigned int i = 0; i <= lhs.degree_; ++i){
        for(unsigned int j = 0; j <= rhs.degree_; ++j){
            product.coefficients_[i + j] += lhs.coefficients_[i] * rhs.coefficients_[j];
        }
    }

    return product;
}

template<typename T, typename TrueType>
void GenericPolynomial<T, TrueType>::pad(unsigned int degree)
{
    if(degree_ >= degree){ // Do nothing if you are already at or above the required degree
        return;
    }

    for(unsigned int i = 0; i < degree - degree_; ++i){ // Add zeroes to the back
        coefficients_.push_back(T());
    }

    degree_ = degree; // Update the degree
}

template<typename T, typename TrueType>
void GenericPolynomial<T, TrueType>::unpad()
{
    int counter = 0;
    for(unsigned int i = degree_; i > 0; --i){
        if(coefficients_[i] == T()){
            ++counter;
        }
        else{
            break;
        }
    }

    degree_ -= counter;
    for(int i = 0; i < counter; ++i){
        coefficients_.pop_back();
    }
}

template<typename T, typename TrueType>
bool GenericPolynomial<T, TrueType>::doubleComp(double a, double b)
{
    return ((a - b) < epsilon_) && ((b - a) < epsilon_);
}

#endif
