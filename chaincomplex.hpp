#ifndef CHAINCOMPLEX_HPP
#define CHAINCOMPLEX_HPP

#include <map>
#include <stdint.h>

#include "homology_z.hpp"

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology > class ChainComplex
{
public:
    ChainComplex();
    MatrixT &operator[] ( int32_t n ) { return differential[n]; }

    Homology homology( int32_t n );
    Homology homology();

private:
    std::map< int32_t, MatrixT& > differential;
};

#include "chaincomplex.ipp"

#endif
