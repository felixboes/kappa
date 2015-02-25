#ifndef DOUBLECOMPLEX_H
#define DOUBLECOMPLEX_H

#include <chrono>
#include <stdint.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_set>

#include <libhomology/homology.hpp>

#include "factorial.hpp"
#include "misc.hpp"
#include "doublecomplexbasis.hpp"
#include "monocomplex.hpp"
#include "operationtester.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"


/**
 * This function is supposed to update the differential according to the contribution
 * of tuple and its boundary to improve readability of the method gen_differential.
 *
 * @warning We don't give a general implementation for this template function since the
 * two specializations we do implement below use specific functions not both of
 * them support.
 *
 * @param parity Parity due to kappa.
 * @param i Number of the differential.
 * @param or_sign Orientation sign.
 */
template <class MatrixType>
void update_differential_doublecomplex(
                         MatrixType &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);
                         
/**
    This DoubleComplex represents the top row.
    One can generate its bases and differentials.
**/
template< class MatrixComplex >
class DoubleComplex
{
public:
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef DoubleComplex< MatrixComplex > ThisType;

    DoubleComplex( const uint32_t genus, const uint32_t num_punctures, SignConvention sgn, const uint32_t number_working_threads, const uint32_t number_remaining_threads);
    
    /**
     *  Recursive function initializing the basis_complex.
    **/
    void gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, Tuple& tuple);
    
    /**
     *  computes the boundary of a given Tuple and saves the result in the differential.
     */ 
    void compute_boundary( Tuple & tuple, const uint32_t p, MatrixType & differential );
    
    /**
     *  Generates the p-th differential.
     */
    void gen_differential( const int32_t p );
    
    /**
     *  @returns a reference to the current differential.
     */
    MatrixType&       get_current_differential();
    
    /**
     *  @returns a reference to the current differential.
     */
    const MatrixType& get_current_differential() const;
    
    /**
     *  @return number of rows of the current differential
    **/
    size_t num_rows() const;
    
    /**
     *  @return number of columns of the current differential
    **/
    size_t num_cols() const;
    
    /**
     *  Access the diagonalizer.
    **/
    DiagonalizerType&       get_diagonalizer();
    
    /**
     *  Access the diagonalizer.
    **/
    const DiagonalizerType& get_diagonalizer() const;
    
    /**
     *  erases the current differential.
     */
    void erase_current_differential();
    
    /**
     *  print a basis to std::out.
     */
    void show_basis( const int32_t p ) const;
//protected:

    uint32_t g;                 ///< genus
    uint32_t m;                 ///< number of punctures
    uint32_t h;                 ///< h = 2*g+m for the parallel case; and h = 2*g+m-1 for the radial case.
    uint32_t num_threads;       ///< number of threads used to construct the differential

    SignConvention sign_conv;   ///< The sign convention.
    MatrixComplex matrix_complex;                         ///< underlying matrix complex of this DoubleComplex
    std::map< int32_t, DoubleComplexBasis > bases;        ///< basis_complex[n] is the n-th DoubleBasis, i.e. the basis of the n-th module of this DoubleComplex. 
};

typedef std::vector<Tuple> DoubleComplexWork;

template< class MatrixComplex >
void doublecomplex_work(DoubleComplex<MatrixComplex> & doublecomplex, DoubleComplexWork & work, const uint32_t p, typename MatrixComplex::MatrixType & differential);

#include "doublecomplex.ipp"

#endif
