#ifndef SAGEMATH_HPP
#define SAGEMATH_HPP

#include <iostream>
#include <sstream>

#include "field_coefficients.hpp"

template< class CoefficientType >
std::string sage_coeff();

template<>
std::string sage_coeff<Q>();

class SagemathInterface
{
public:
    SagemathInterface();
    ~SagemathInterface();

    void test() const;

    template< class CoefficientType >
    void create_matrix( size_t rows, size_t cols );

    template< class MatrixType >
    void update_row( const int32_t row, const MatrixType& vect );

    void compute_rank();

    FILE* sagemath_pipe;
};

#include "sagemath.ipp"

#endif // SAGEMATH_HPP
