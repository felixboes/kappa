#include "matrix.hpp"

template< class CoefficientT >
Matrix<CoefficientT>::Matrix() : data(), num_cols(0), num_rows(0)
{
}

template< class CoefficientT >
Matrix<CoefficientT>::Matrix( size_t number_rows, size_t number_cols ) : data( number_rows*number_cols, CoefficientT(0) ), num_rows(number_rows), num_cols(number_cols)
{
}

template< class CoefficientT >
void Matrix<CoefficientT>::row_operation( size_t row_1, size_t row_2, size_t col )
{
    CoefficientT lambda( -at(row_1,col) / at(row_2, col) );
    for( size_t j = col; j < num_cols; ++j )
    {
        CoefficientT & a = this->operator()( row_2, j );
        a = lambda * a + at( row_1, j );
    }
}

template< class CoefficientT >
 CoefficientT & Matrix<CoefficientT>::operator()( size_t i, size_t j )
 {
     return data[ i*num_cols + j ];
 }

template< class CoefficientT >
const CoefficientT & Matrix<CoefficientT>::at(size_t i, size_t j) const
{
    return data.at( i*num_cols + j );
}

template< class CoefficientT >
size_t Matrix<CoefficientT>::size1() const
{
    return num_rows;
}

template< class CoefficientT >
size_t Matrix<CoefficientT>::size2() const
{
    return num_cols;
}
