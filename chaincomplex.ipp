#include "chaincomplex.hpp"

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology > ChainComplex::ChainComplex()
{


}

template< class Coefficient, class MatrixT, class DiagonalizerT, class Homology> Homology ChainComplex::homology(int32_t n)
{
    DiagonalizerT Diag( module[n+1], module[n] );
    return Diag.homology();
}
