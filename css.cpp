#include "css.hpp"

uint32_t CSSBasis :: add_basis_element (Tuple& t)
{
    uint32_t num_clusters = t.num_cycles();
    LBasisType l_basis = basis[num_clusters];
    t.id = l_basis.size();
    l_basis.insert(t);
    
    return t.id;
}

uint64_t CSSBasis :: size( uint32_t l ) const
{
    return (basis.count(l) != 0 ? basis.at(l).size() : 0 );
}

uint64_t CSSBasis :: total_size() const
{
    uint64_t size(0);
    for( auto& it : basis )
    {
        size += it.second.size();
    }
    return size;
}

int64_t CSSBasis :: id_of(Tuple &t) const
{
    auto it = basis.find(t);
    if( it == basis.end() )
    {
        return -1;
    }
    else
    {
        return it->id;
    }
}

std::ostream& operator<< ( std::ostream& os, const CSSBasis& cb )
{
    for( auto& it : cb.basis )
    {
        os << "Cluster of size " << it->first << std::endl;
        for( auto& b : it.second )
        {
            os << b.id << ": " << b << std::endl;
        }
    }
    return os;
}
