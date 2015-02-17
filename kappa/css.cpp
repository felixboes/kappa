#include "css.hpp"

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
