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


#ifndef CSS_H
#define CSS_H

#include <algorithm>
#include <chrono>
#include <stdint.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <unordered_set>
#include <libhomology/homology.hpp>

#include "cssbasis.hpp"
#include "factorial.hpp"
#include "monocomplex.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"


/**
This ClusterSpectralSequence represents E_1 term of the cluster spectral sequence associated with the Ehrenfried complex 
which is generated by the monotone tuples of transpositions. 
**/

template< class MatrixComplex >
class ClusterSpectralSequence
{
public:
    typedef typename MatrixComplex::CoefficientType CoefficientType;
    typedef typename MatrixComplex::MatrixType MatrixType;
    typedef typename MatrixComplex::DiagonalizerType DiagonalizerType;
    typedef typename MatrixComplex::HomologyType HomologyType;
    typedef std::map< int32_t, HomologyType > CSSHomologyType;
    
    ClusterSpectralSequence(
            const uint32_t          genus,
            const uint32_t          num_punctures,
            const SignConvention    sgn,
            const uint32_t          number_working_threads,
            const uint32_t          number_remaining_threads );
    /** Recursive function initializing the basis_complex.
        In the call of gen_bases with the parameters s, p and tuple, we assume that the first s transpositions
        containing symbols 1, ..., p are fixed and append all possible transpositions at position s+1, applying
        the function recursively in an appropriate way.
        If s == h, we don't append another transposition since we have completed a possible basis element. 
        We check whether its number of cycles is appropriate for it to be a basis element, and if this is
        the case, we add it to the basis in degree p. Thereby, basis elements are sorted according to the number
        of clusters.
    **/
    void gen_bases( const uint32_t l, const uint32_t p, const uint32_t start_symbol, Tuple& tuple);
    void gen_d0( int32_t p, int32_t l );
    void gen_d0_boundary(const Tuple & tuple,
                         const int32_t p,
                         const int32_t l,
                         typename MatrixComplex::MatrixType & differential);
    
    void gen_d1_stage_1( const int32_t p, const int32_t l );
    MatrixType gen_d1_row( const int32_t, const int32_t l, const Tuple& basis_element );
    void gen_d1_apply_operations( MatrixType& row );
    void prepare_d1_diag();
    void erase_d0();
    void erase_d1();
    
    void show_basis( const int32_t p ) const;         ///< print a basis to std::out
//protected:

    uint32_t g;     ///< genus
    uint32_t m;     ///< number of punctures
    uint32_t h;     ///< h = 2*g+m
    uint32_t num_threads;  ///< number of threads used to construct the differential
    
    SignConvention sign_conv;                    ///< The sign convention.
    MatrixComplex diff_complex;                  ///< Due to RAM limitations, we are working with at most two matrices at a time. Therefore we do not model the whole spectral sequence.
    std::map< int32_t, CSSBasis > basis_complex; ///< basis_complex[n] is the n-th CSSBasis.
    
};

typedef std::vector<Tuple> CSSWork;

template< class MatrixComplex >
void css_work_0(ClusterSpectralSequence<MatrixComplex> & css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              typename MatrixComplex::MatrixType & differential
              );

template< class MatrixComplex >
void css_work_1(ClusterSpectralSequence<MatrixComplex> & css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              typename MatrixComplex::MatrixType & differential,
              const size_t num_cols,
              const std::vector< size_t >& offset
              );

// Template specializations for bool-matrices.
template<>
void ClusterSpectralSequence< ChainComplexBoolCSS >::gen_d1_apply_operations( MatrixType& row );

template<>
void css_work_1(ClusterSpectralSequence< ChainComplexBoolCSS >& css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              ChainComplexBoolCSS::MatrixType& differential,
              const size_t num_cols,
              const std::vector< size_t >& offset
              );

#endif // CSS_HPP
