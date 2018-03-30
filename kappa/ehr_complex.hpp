// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#ifndef EHR_COMPLEX_H
#define EHR_COMPLEX_H

#include <chrono>
#include <stdint.h>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_set>
#include <sstream>

#include <libhomology/homology.hpp>

#include "factorial.hpp"
#include "misc.hpp"
#include "ehr_basis.hpp"
#include "sessionconfig.hpp"
#include "sym_grp_tuple.hpp"
#include "alt_grp_tuple.hpp"


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
void update_differential(MatrixType &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);

/**
 * This template specialization updates the differential of type MatrixBool according to the contribution
 * of tuple and its boundary.
 */
template<>
void update_differential(MatrixBool &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);

template<>
void update_differential(MatrixBoolCSS &        differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          parity,
                         const int8_t           i,
                         const int8_t           or_sign,
                         const SignConvention & sign_conv);
                         
/**
    This EhrComplex represents a chain complex which is generated by the fully unstable tuples of transpositions.
    The EhrComplex can either consist of parallel or of radial cells, which is marked by a flag.
    One can generate its bases and differentials.
**/
template< class MatrixComplex, class TupleT >
class EhrComplex
{
public:
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef TupleT TupleType;
    typedef typename TupleT::HashTuple HashTupleType;
    typedef typename TupleT::EhrBasesGenerator EhrBasesGeneratorType;
    typedef EhrComplex< MatrixComplex, TupleT > ThisType;

    EhrComplex(
            const uint32_t          genus,
            const uint32_t          num_punctures,
            const SignConvention    sgn,
            const uint32_t          number_working_threads,
            const uint32_t          number_remaining_threads );
    
    /**
     *  computes the boundary of a given Tuple and saves the result in the differential.
     */ 
    void compute_boundary( TupleT & tuple, const uint32_t p, MatrixType & differential);
    
    /**
     *  Generates the p-th differential.
     */
    void gen_differential( const int32_t p );
    
    /**
     *  Apply base changes.
    **/
    void apply_base_changes();
    
    /**
     *  Compute Homchain.
    **/
    void homchain(int32_t p, bool homology = false, int32_t maxdimension=-1);

    /**
     *  Diagoanlize current differential.
     */
    HomologyType diagonalize_current_differential( const int32_t p, uint32_t max_rank = 0, const bool print_duration = true );
    
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
     *  Compute the kernel at the \f$n\f$-th spot and the torsion at the \f$(n-1)\f$-th spot.
    **/
    HomologyType compute_current_kernel_and_torsion( const int32_t n );
    
    /**
     *  print a basis to std::out.
     */
    void show_basis( const int32_t p ) const;
//protected:

    uint32_t g;                 ///< genus
    uint32_t m;                 ///< number of punctures
    uint32_t h;                 ///< h = 2*g+m for the parallel case; and h = 2*g+m-1 for the radial case.
    uint32_t num_threads;       ///< number of threads used to construct the differential

    SignConvention sign_conv;                    ///< The sign convention.
    MatrixComplex diff_complex;                  ///< Due to RAM limitations, we are working with at most two matrices at a time. Therefore we do not model the whole spectral sequence.
    std::map< int32_t, EhrBasis<TupleT> > basis_complex; ///< basis_complex[n] is the n-th EhrBasis.
    std::ofstream homchain_homology_file;
    std::ofstream homchain_cohomology_file;
};

template<class TupleT>
using EhrComplexWork = std::vector<TupleT>;

template< class MatrixComplex, class TupleT >
void ehr_complex_work(EhrComplex<MatrixComplex, TupleT> &ehrcomplex, EhrComplexWork<TupleT> &work, const uint32_t p,
                      typename MatrixComplex::MatrixType &differential);

/**
 * @brief Formula for the sign appearing in d_hor.
 */
int32_t sign(const int32_t          parity,
             const int8_t           i,
             const int8_t           or_sign,
             const SignConvention & sign_conv );

#endif // EHR_COMPLEX_H
