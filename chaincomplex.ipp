#include "chaincomplex.hpp"

template< class Coefficient, class MatrixT, class DiagonalizerT > ChainComplex::ChainComplex()
{


}

template< class Coefficient, class MatrixT, class DiagonalizerT > HomologyModule ChainComplex::homology(int32_t n)
{
    DiagonalizerT Diag( module[n+1], module[n] );
    return Diag.homology();
}
