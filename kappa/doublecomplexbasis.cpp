#include "doublecomplexbasis.hpp"

DoubleComplexBasis::DoubleComplexBasis() : basis_h(), basis_h_1()
{
}
 
uint32_t DoubleComplexBasis :: add_basis_element ( HighCell t )
{
    uint32_t id = 0;
    if( t.is_redundant() )
    {
        id = basis_h_1.size();
        t.id = id;
        basis_h_1.insert(std::move(t));
    }
    else
    {
        id = basis_h.size();
        t.id = id;
        basis_h.insert(std::move(t));
    }
    
    return id;
}

uint64_t DoubleComplexBasis :: size_h() const
{
    return basis_h.size();
}

uint64_t DoubleComplexBasis :: size_h_1() const
{
    return basis_h_1.size();
}

int64_t DoubleComplexBasis :: id_of(const HighCell &t) const
{
    if( t.is_redundant() )
    {
        auto it = basis_h_1.find(t);
        if( it == basis_h_1.end() )
        {
            return -1;
        }
        else
        {
            return it->id;
        }
    }
    else
    {
        auto it = basis_h.find(t);
        if( it == basis_h.end() )
        {
            return -1;
        }
        else
        {
            return it->id;
        }
    }
}

std::ostream& operator<< ( std::ostream& os, const DoubleComplexBasis& dcb )
{
    os << "essential and collapsable of horizontal: " << std::endl;
    for( const auto& it : dcb.basis_h )
    {
        os << it.id << ": " << it << std::endl;
    }
    os << "redundant: " << std::endl;
    for( const auto& it : dcb.basis_h_1 )
    {
        os << it.id << ": " << it << std::endl;
    }
    
    return os;
}
