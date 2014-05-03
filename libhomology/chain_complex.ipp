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
const MatrixT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::at (const int32_t& n ) const
{
    return differential.at(n);
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
size_t ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::count (const int32_t& n ) const
{
    return differential.count(n);
}

template < class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
void ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::erase (const int32_t &n)
{
    if( differential.count(n) != 0 )
    {
        // Delete matrix to save space.
        // Quote from the boost::ublas documentation http://www.boost.org/doc/libs/1_49_0/libs/numeric/ublas/doc/matrix.htm
        //
        // void resize (size_type size1, size_type size2, bool preserve = true)
        // Reallocates a matrix to hold size1 rows of size2 elements. The existing elements of the matrix are preseved when specified.
        differential[n].resize(0,0,false);
        differential.erase(n);
    }
}

template < class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
void ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::erase_all ()
{
    for( auto& diff_it: differential )
    {
        diff_it.second.resize(0,0,false);
    }
    differential.clear();
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
CoefficientT &ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::operator() ( int32_t n, uint32_t row, uint32_t col )
{
    return differential[n](row, col);
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
bool ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::exists_differential( const int32_t& n ) const
{
    // The n-th module exists iff there exists a matrix leaving it.
    return differential.count(n) > 0;
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
HomologyT ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::homology( int32_t n, uint32_t number_threads )
{
    atomic_uint current_rank;
    return homology(n, current_rank, number_threads);
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
            typedef typename HomologyT::KernT KT;
            KT k = KT();
            homol.set_kern( n, k );
            std::cerr << "Error: Homology in " << n << " unknown." << std::endl;
        }
    }
    // diff_n has no rows.
    else if ( differential[n].size1() == 0 )
    {
        // The dimension equals the number of columns.
        homol.set_kern( n, differential[n].size2() );
        // There is no torsion one dimension below.
        typedef typename HomologyT::TorsT TT;
        TT t = TT();
        homol.set_tors( n-1, t);
    }
    else
    {
        // Diagonalize.
        DiagonalizerT diago;
        diago( differential[n], current_rank, number_threads );
        homol.set_kern( n, diago.kern() );
        homol.set_tors( n-1, diago.tors() );
    }
    
    return homol;    
}

template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
HomologyT ChainComplex< CoefficientT, MatrixT, DiagonalizerT, HomologyT >::homology( int32_t n, atomic_uint & current_rank, uint32_t number_threads )
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
            typedef typename HomologyT::KernT KT;
            KT k = KT();
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
        diago( differential[n], current_rank, number_threads );
        homol.set_kern( n, diago.kern() );
    }
    // diff_{n+1} is 0
    if( differential.count(n+1) == 0 || differential[n+1].size1() == 0 )
    {
        // There is no torsion.
        typedef typename HomologyT::TorsT TT;
        TT t = TT();
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
        diago( differential[n+1], current_rank, number_threads );
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
    for( auto& it : differential )
    {
        int32_t n = it.first;
        MatrixT diff = it.second;
        // if diff_n is 0, everything is kernel and nothing is torsion
        if( diff.size1() == 0 )
        {
            homol.set_kern( n, diff.size2() );
            typedef typename HomologyT::TorsT TT;
            TT t = TT();
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
