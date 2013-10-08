#ifndef CHAIN_COMPLEX_HPP
#define CHAIN_COMPLEX_HPP

#include <cinttypes>
#include <iostream>
#include <map>

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology >
class ChainComplex
{
public:
    ChainComplex();
    
    /**
     *  Access the \f$n\f$-th differential.
     */
    MatrixT &operator[] ( int32_t n );
    
    /**
     *  Access the coefficient of the \f$n\f$-th differential at the position (row, col).
     */
    Coefficient &operator() ( int32_t n, uint32_t row, uint32_t col );
    
    /**
     *  Compute the homology at the \f$n\f$-th spot.
     */
    Homology homology( int32_t n );
    
    /**
     *  Compute all the homology.
     */
    Homology homology();

    /**
     *  Checks if the \f$n\f$-th modules exists.
     */
    bool exists_module( int32_t n );
    
private:
    std::map< int32_t, MatrixT > differential;
};

#include "chain_complex.ipp"

#endif
