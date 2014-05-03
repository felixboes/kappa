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
  *  Therefore we imagine a chaincomplex as a list of matrices.
  *  In order to avoid a high memory usage and reallocations,
  *  we only store a single differnential at each time.
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
    HomologyT compute_kernel_and_torsion( int32_t n, uint32_t number_threads = 1 );
    
    /**
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
     */
    HomologyT compute_kernel_and_torsion( int32_t n, atomic_uint & current_rank, uint32_t number_threads = 1 );
    
private:
    MatrixT current_differential; ///< Realizes the transpose of a single differential.

    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & current_differential;
    }
};

#include "chain_complex.ipp"

#endif
