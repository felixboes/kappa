#include "chain_complex.hpp"

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::ChainComplex()
{
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
MatrixT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::get_current_differential( int32_t n )
{
    return current_differential;
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
const MatrixT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::get_current_differential( int32_t n ) const
{
    return current_differential;
}

template < class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
void ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::erase ()
{
    current_differential.erase();
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
size_t ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::num_cols() const
{
    return current_differential.size2();
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
size_t ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::num_rows() const
{
    return current_differential.size1();
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
CoefficientT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::operator() ( uint32_t row, uint32_t col )
{
    return current_differential(row, col);
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
HomologyT ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::compute_kernel_and_torsion( int32_t n, uint32_t number_threads )
{
    atomic_uint current_rank(0);
    return compute_kernel_and_torsion(n, current_rank, number_threads);
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
HomologyT ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::compute_kernel_and_torsion( int32_t n, atomic_uint & current_rank, uint32_t number_threads )
{
 {
    HomologyT homol;
    // Diagonalize.
    DiagonalizerT diago;
    diago( differential[n], current_rank, number_threads );
    homol.set_kern( n, diago.kern() );
    homol.set_tors( n-1, diago.tors() );
    
    return homol;    
}

