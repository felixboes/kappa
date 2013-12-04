#include "chain_complex.hpp"

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::ChainComplex()
{
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
MatrixT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::operator[] ( int32_t n )
{
    return differential[n];
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
const MatrixT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::at (const int32_t& n )
{
    return differential.at(n);
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
CoefficientT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::operator() ( int32_t n, uint32_t row, uint32_t col )
{
    return differential[n](row, col);
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
bool ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::exists_differential( int32_t n )
{
    // The n-th module exists iff there exists a matrix leaving it.
    return differential.count(n) > 0;
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
HomologyT ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::homology( int32_t n )
{
    HomologyT homol;
    // Case by case anaylsis:
    
    // diff_n is the zero matrix if either it is not stored or it has no rows.
    // In this case, the defect is the whole module.
    
    // diff_n is not stored.
    if( differential.count(n) == 0 )
    {
        // The dimension of the n-th module is known by the (n+1)-th differential.
        if( differential.count(n+1) )
        {
            homol.set_kern( n, differential[n+1].size1() );
        }
        else // The dimension of the n-th module is unkown.
        {
            typename HomologyT::KernT k;
            homol.set_kern( n, k );
            std::cerr << "Error: Homology in " << n << " unknown." << std::endl;
        }
    }
    // diff_n has no rows.
    else if ( differential[n].size1() == 0 )
    {
        // The dimension equals the number of columns.
        homol.set_kern( n, differential[n].size2() );
    }
    else
    {
        // Diagonalize.
        DiagonalizerT diago;
        diago( differential[n] );
        homol.set_kern( n, diago.kern() );
    }
    // diff_{n+1} is 0
    if( differential.count(n+1) == 0 || differential[n+1].size1() == 0 )
    {
        // Tthere is no torsion.
        typename HomologyT::TorsT t;
        homol.set_tors( n, t );
    }
    else 
    {
        // Assert that the matrices have the right number of colums and rows.
        if( !(
                differential.count(n) == 0 || differential[n].size1() == 0 || // Zero matrix
                differential[n].size2() == differential[n+1].size1()
                ) )
        {
            std::cerr << "Error: Column-Rows exception at position " << n << std::endl;
        }
        // Diagonalize.
        DiagonalizerT diago;
        diago( differential[n+1] );
        homol.set_tors( n, diago.tors() );
    }
    
    return homol;
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
HomologyT ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::homology()
{
    HomologyT homol;
    
    // Iterate through all the existing differentials. For each differential diff_n, we
    // compute the kernel (stored in the n-th homology module) and the image (stored
    // as torsion in the (n-1)-st homology module).
    for( auto it = differential.begin(); it != differential.end(); ++it )
    {
        int32_t n = it->first;
        MatrixT diff = it->second;
        // if diff_n is 0, everything is kernel and nothing is torsion
        if( diff.size1() == 0 )
        {
            homol.set_kern( n, diff.size2() );
            typename HomologyT::TorsT t(0);
            homol.set_tors( n-1, t );
        }
        else
        {
            // diagonalize
            DiagonalizerT diago;
            diago(diff);
            homol.set_kern( n, diago.kern() );
            homol.set_tors( n-1, diago.tors() );
        }
    }
    
    return homol;
}

