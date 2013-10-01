#ifndef CHAINCOMPLEX_HPP
#define CHAINCOMPLEX_HPP

#include <map>
#include <stdint.h>

#include "homology.hpp"

template< class Coefficient, class MatrixT, class DiagonalizerT > class ChainComplex
{
public:
    ChainComplex();
    MatrixT &operator[] ( int32_t n ) { return differential[n]; }

    HomologyModule< Coefficient > homology( int32_t n );
    HomologyComplex< Coefficient > homology();

private:
    std::map< int32_t, MatrixT& > differential;
};

#include "chaincomplex.ipp"

#endif
