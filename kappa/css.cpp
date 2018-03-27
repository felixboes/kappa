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


#include "css.hpp"
#include "css_impl.ipp"
#include "ehr_complex_impl.ipp"

template<>
void ClusterSpectralSequence< ChainComplexBoolCSS >::gen_d1_apply_operations( MatrixType& row )
{
    MatrixType& differential = diff_complex.get_current_differential();
    const typename MatrixType::DiagonalType& diagonal = differential.diagonal;
    
    for( const auto& diag_pos : diagonal )
    {
        const auto& diag_row = diag_pos.first;
        const auto& diag_col = diag_pos.second;
        
        if( row.at( 0, diag_col ) != 0 )
        {
            row.main_row(0) ^= differential.main_row_at(diag_row);
        }
    }
}

template<>
void css_work_1(ClusterSpectralSequence< ChainComplexBoolCSS >& css,
              CSSWork & work,
              const int32_t p,
              const int32_t l,
              ChainComplexBoolCSS::MatrixType& differential,
              const size_t num_cols,
              const std::vector< size_t >& offset)
{
    for ( auto it : work)
    {
        MatrixBoolCSS single_row = css.gen_d1_row( p, l, it );
        css.gen_d1_apply_operations( single_row );
        for( size_t j = 0; j < num_cols; ++j )
        {
            differential.sec_set( it.id, j, single_row.at( 0, j + offset[j] ) );
        }
    }
}

#define force_template_instantiation( MatrixComplex, TupleType )\
    template class ClusterSpectralSequence<MatrixComplex>;\
    template void update_differential(MatrixComplex &differential, const size_t row, const size_t column, const int32_t, const int8_t, const int8_t, const SignConvention &);\
    template void ehr_complex_work(EhrComplex<MatrixComplex, TupleType> & ehrcomplex, EhrComplexWork<TupleType> & work, const uint32_t p, MatrixComplex::MatrixType & differential);\
    template void css_work_0(ClusterSpectralSequence< MatrixComplex > & css, CSSWork & work, const int32_t p, const int32_t l, typename MatrixComplex::MatrixType & differential );\
    template void css_work_1(ClusterSpectralSequence< MatrixComplex > & css, CSSWork & work, const int32_t p, const int32_t l, typename MatrixComplex::MatrixType & differential, const size_t num_cols, const std::vector< size_t >& offset );

force_template_instantiation( ChainComplexQCSS, SymGrpTuple )
force_template_instantiation( ChainComplexZmCSS, SymGrpTuple )

#undef force_template_instantiation
