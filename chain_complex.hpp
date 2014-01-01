#ifndef CHAIN_COMPLEX_HPP
#define CHAIN_COMPLEX_HPP

#include <cinttypes>
#include <iostream>
#include <map>

#include <boost/serialization/access.hpp>

#include "parallelization.hpp"

/**
 *  We assume, that \f$ 0 \times n \f$- and \f$ n \times 0 \f$-matrices require negligible RAM.
 *  Using standard tools, our tests where successfull for the types MatrixZm and MatrixQ.
 */
template< class CoefficientT, class MatrixT, class DiagonalizerT, class HomologyT >
class ChainComplex
{
public:
    typedef ChainComplex<CoefficientT, MatrixT, DiagonalizerT, HomologyT> SelfType;
    typedef CoefficientT CoefficientType;
    typedef MatrixT MatrixType;
    typedef DiagonalizerT DiagonalizerType;
    typedef HomologyT HomologyType;
    
    ChainComplex();
    
    
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
    HomologyT homology( int32_t n );
    
    /**
     *  Compute the homology at the \f$n\f$-th spot.
     */
    HomologyT homology( int32_t n, atomic_uint & current_rank );
    
    /**
     *  Compute all the homology.
     */
    HomologyT homology();

    /**
     *  Checks if the \f$n\f$-th differential exists.
     */
    bool exists_differential( int32_t n );
    
private:
    std::map< int32_t, MatrixT > differential;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
        ar & differential;
    }
};

#include "chain_complex.ipp"

#endif
