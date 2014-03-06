#ifndef CHAIN_COMPLEX_HPP
#define CHAIN_COMPLEX_HPP

// Description:
//
// This header defines a very generic chain complex type.
// Thus the implementation may look a bit messy at first glance.
// Our chain complex depents on coefficients, matrix types, a function object that diagonalizes such matrices and a homology type.

#include <cinttypes>
#include <iostream>
#include <map>

#include <boost/serialization/access.hpp>

#include "parallelization.hpp"

/**
 *  In our applications we do not always need explicit basis.
 *  Therefore we realized a chaincomplex as a list of matrices.
 *
 *  We assume, that \f$ 0 \times n \f$- and \f$ n \times 0 \f$-matrices require negligible RAM.
 */
template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
class ChainComplex
{
public:
    typedef ChainComplex<CoefficientT, MatrixT, DiagonalizerT, HomologyT> SelfType; ///< We use this typedef to grant access this type from other classes.
    typedef CoefficientT CoefficientType;   ///< We use this typedef to grant access this type from other classes.
    typedef MatrixT MatrixType;             ///< We use this typedef to grant access this type from other classes.
    typedef DiagonalizerT DiagonalizerType; ///< We use this typedef to grant access this type from other classes.
    typedef HomologyT HomologyType;         ///< We use this typedef to grant access this type from other classes.
    
    ChainComplex(); ///< Constructs an empty chain complex.
    
    /**
     *  Access the \f$n\f$-th differential.
     */
    MatrixT &operator[] ( int32_t n );
    
    /**
     *  Access the \f$n\f$-th differential.
     */
    const MatrixT &at ( const int32_t& n ) const;
    
    /**
     *  @returns 0 iff there is no differential stored and 1 else.
     */
    size_t count( const int32_t& n ) const;
    
    /**
     *  Erases the n-th differential of the complex.
     */
    void erase( const int32_t& n );
    
    /**
     *  Access the coefficient of the \f$n\f$-th differential at the position (row, col).
     */
    CoefficientT &operator() ( int32_t n, uint32_t row, uint32_t col );
    
    /**
     *  Compute the homology at the \f$n\f$-th spot.
     */
    HomologyT homology( int32_t n, uint32_t number_threads = 0 );
    
    /**
     *  Compute the homology at the \f$n\f$-th spot.
     */
    HomologyT homology( int32_t n, atomic_uint & current_rank, uint32_t number_threads = 0 );

    /**
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
     */
    HomologyT compute_kernel_and_torsion( int32_t n, uint32_t number_threads=0 );
    
    /**
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
     */
    HomologyT compute_kernel_and_torsion( int32_t n, atomic_uint & current_rank, uint32_t number_threads=0 );
    
    /**
     *  Compute all the homology.
     */
    HomologyT homology();

    /**
     *  Checks if the \f$n\f$-th differential exists.
     */
    bool exists_differential( const int32_t& n ) const;
    
private:
    std::map< int32_t, MatrixT > differential;  ///< Realizes the data.
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & differential;
    }
};

#include "chain_complex.ipp"

#endif
