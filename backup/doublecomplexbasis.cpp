#include "doublecomplexbasis.hpp"

DoubleComplexBasis::DoubleComplexBasis() : basis_red(), basis_col(), basis_ess()
{
}
 
void DoubleComplexBasis :: add_basis_element ( HighCell t )
{
    if( t.is_redundant() )
    {
        basis_red.insert(std::move(t));
    }
    else if( t.monotone() )
    {
        basis_ess.insert(std::move(t));
    }
    else
    {
        basis_col.insert(std::move(t));
    }
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

void DoubleComplexBasis::generate_indices()
{
    size_t id = 0;
    for( auto& cell : basis_red )
    {
        cell.id = id;
        ++id;
    }
    id = 0;
    for( auto& cell : basis_col )
    {
        cell.id = id;
        ++id;
    }
    id = 0;
    for( auto& cell : basis_ess )
    {
        cell.id = id;
        ++id;
    }
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
