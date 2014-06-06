#include "css.hpp"

int32_t CSSBasis :: add_basis_element ( Tuple& t )
{
    uint32_t num_clusters = t.num_cluster();
    LBasisType& l_basis = basis[num_clusters];
    t.id = l_basis.size();
    l_basis.insert(t);
    
    return t.id;
}

int32_t CSSBasis :: size( const int32_t l ) const
{
    return (basis.count(l) != 0 ? basis.at(l).size() : 0 );
}

int32_t CSSBasis :: total_size() const
{
    int32_t size(0);
    for( auto& it : basis )
    {
        size += it.second.size();
    }
    return size;
}

int32_t CSSBasis :: id_of( Tuple& t ) const
{
    for( auto& l_basis_it : basis )
    {
        auto& l_basis = l_basis_it.second;
        auto it = l_basis.find(t);
        if( it != l_basis.end() )
        {
            return it->id;
        }
    }
    return -1;
}

int32_t CSSBasis :: total_id_of( Tuple& t ) const
{
    int64_t basis_offset(0);
    
    for( auto& l_basis_it : basis )
    {
        auto& l_basis = l_basis_it.second;
        auto it = l_basis.find(t);
        if( it != l_basis.end() )
        {
            return basis_offset + it->id;
        }
        basis_offset += this->size(l_basis_it.first);
    }
    return -1;
}

std::ostream& operator<< ( std::ostream& os, const CSSBasis& cb )
{
    for( auto& it : cb.basis )
    {
        os << "Cluster of size " << it.first << std::endl;
        for( auto& b : it.second )
        {
            os << b.id << ": " << b << std::endl;
        }
    }
    return os;
}

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
