#include "doublecomplexbasis.hpp"

DoubleComplexBasis::DoubleComplexBasis() : basis_red(), basis_col(), basis_ess()
{
}
 
uint32_t DoubleComplexBasis :: add_basis_element ( HighCell t )
{
    uint32_t id = 0;
    if( t.is_redundant() )
    {
        id = basis_red.size();
        t.id = id;
        basis_red.insert(std::move(t));
    }
    else if( t.monotone() )
    {
        id = basis_ess.size();
        t.id = id;
        basis_ess.insert(std::move(t));
    }
    else
    {
        id = basis_col.size();
        t.id = id;
        basis_col.insert(std::move(t));
    }
    
    return id;
}

uint64_t DoubleComplexBasis :: size_red() const
{
    return basis_red.size();
}

uint64_t DoubleComplexBasis :: size_col() const
{
    return basis_col.size();
}

uint64_t DoubleComplexBasis :: size_ess() const
{
    return basis_ess.size();
}

int64_t DoubleComplexBasis :: id_of(const HighCell &t) const
{
    if( t.is_redundant() )
    {
        auto it = basis_red.find(t);
        if( it == basis_red.end() )
        {
            return -1;
        }
        else
        {
            return it->id;
        }
    }
    else if( t.monotone() )
    {
        auto it = basis_ess.find(t);
        if( it == basis_ess.end() )
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
        auto it = basis_col.find(t);
        if( it == basis_col.end() )
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
    if( dcb.size_col() > 0 )
    {
        os << "Collapsible cells: " << std::endl;
        for( const auto& it : dcb.basis_col )
        {
            os << it.id << ": " << it << std::endl;
        }
    }
    else
    {
        os << "There are no collapsible cells." << std::endl;
    }

    if( dcb.size_ess() > 0 )
    {
        os << "Essential cells: " << std::endl;
        for( const auto& it : dcb.basis_ess )
        {
            os << it.id << ": " << it << std::endl;
        }
    }
    else
    {
        os << "There are no essential cells." << std::endl;
    }
    
    if( dcb.size_red() > 0 )
    {
        os << "Redundant cells: " << std::endl;
        for( const auto& it : dcb.basis_red )
        {
            os << it.id << ": " << it << std::endl;
        }
    }
    else
    {
        os << "There are no redundant cells." << std::endl;
    }
    
    return os;
}
