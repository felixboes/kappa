#include "chain_complex.hpp"

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology >
ChainComplex< Coefficient, MatrixT, DiagonalizerT, Homology >::ChainComplex()
{


}

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology >
MatrixT &ChainComplex< Coefficient, MatrixT, DiagonalizerT, Homology >::operator[] ( int32_t n )
{
    return differential[n];
}

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology >
Coefficient &ChainComplex< Coefficient, MatrixT, DiagonalizerT, Homology >::operator() ( int32_t n, uint32_t row, uint32_t col )
{
    return differential[n](row, col);
}

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology >
bool ChainComplex< Coefficient, MatrixT, DiagonalizerT, Homology >::exists_module( int32_t n )
{
    // The n-th module exists iff there exists a matrix M leaving it.
    return differential.count(n) > 0;
}

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology >
Homology ChainComplex< Coefficient, MatrixT, DiagonalizerT, Homology >::homology( int32_t n )
{
    if( exists_module(n+1) && exists_module(n) )
    {
        DiagonalizerT diag( differential[n+1], differential[n] );
        return Homology( n, diag.kern(), diag.torsion() );
    }
    else
    {
        std::cerr << "Todo" << std::endl;
    }
}
