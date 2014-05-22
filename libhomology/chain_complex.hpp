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
  *  Therefore we imagine a chaincomplex as a list of matrices.
  *  There are two possibilities to use this class:
  *
  *  1) If one is interested in computing the homology of a chain complex
  *     with huge differentials, one can use the data member current_differential
  *     to store only one differential of the chain complex at each time.
  *     This avoids a high memory usage and reallocations.
  *     The methods in the first sections are provided for this approach.
  *  2) One can store the whole chain complex in the map differential, where
  *     for each reasonable p, the pth differential is stored. This approach is useful
  *     for small chain complex, or if you want to re-use the base change of the nth
  *     differential for diagonalizing the (n+1)th differential.
  *     The second section provides methods for this approach.
  *
  *  We assume that \f$ 0 \times n \f$- and \f$ n \times 0 \f$-matrices require negligible RAM.
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
    
    ChainComplex(bool matrices_are_transposed = false) : transp(matrices_are_transposed) {} ///< Constructs an empty chain complex.

    // Methods for the usage of current_differential.

     /**
      *  Access the transpose of the current differential.
      */
     MatrixT &       get_current_differential();
     const MatrixT & get_current_differential() const;
     
     /**
      *  Erases the current differential of the complex.
      */
     void erase();
    
     /**
      *  Access the coefficient of the current differential at the position (row, col).
      */
     CoefficientT &operator() ( uint32_t row, uint32_t col );
     
     /**
       * @return number of rows resp. columns of the current differential
       */
     size_t num_rows() const;
     size_t num_cols() const;

    /**
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
     */
    HomologyT compute_current_kernel_and_torsion( int32_t n, uint32_t number_threads = 1 );
    
    /**
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
     */
    HomologyT compute_current_kernel_and_torsion( int32_t n, atomic_uint & current_rank, uint32_t number_threads = 1 );

/////////////////////////////////////////////////////////////////////////////////////////////////////

    // Methods for the use of the map of differentials.

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
     *  Erases all differentials of the complex.
     */
    void erase_all();

    /**
     *  Access the coefficient of the \f$n\f$-th differential at the position (row, col).
     */
    CoefficientT &operator() ( int32_t n, uint32_t row, uint32_t col );

    /**
     *  Compute the homology at the \f$n\f$-th spot.
     */
    HomologyT homology( int32_t n, uint32_t number_threads = 1 );

    /**
    *  Compute the homology at the \f$n\f$-th spot.
    */
    HomologyT homology( int32_t n, atomic_uint & current_rank, uint32_t number_threads = 0 );

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
    bool transp;
    MatrixT current_differential; ///< Realizes the transpose of a single differential.

    std::map< int32_t, MatrixT > differential;  ///< Realizes all differentials.

    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & current_differential;
        ar & differential;
    }
};

#include "chain_complex.ipp"

#endif
