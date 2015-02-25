#include "doublecomplexbasis.hpp"

DoubleComplexBasis::DoubleComplexBasis() : basis()
{
}
 
uint32_t DoubleComplexBasis :: add_basis_element ( Tuple t)
{
    uint32_t id = basis.size();
    t.id = id;
    basis.insert(std::move(t));
    
    return id;
}

uint64_t DoubleComplexBasis :: size() const
{
    return basis.size();
}

int64_t DoubleComplexBasis :: id_of(const Tuple &t) const
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

std::ostream& operator<< ( std::ostream& os, const DoubleComplexBasis& dcb )
{
    for( const auto& it : dcb.basis )
    {
        os << it.id << ": " << it << std::endl;
    }
    return os;
}
