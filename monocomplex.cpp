#include "monocomplex.hpp"

/*
 *
 *   MonoBasis
 *
 */
uint32_t MonoBasis :: add_basis_element (Tuple& t)
{
    t.id = basis.size();
    basis.insert(t);
    
    return t.id;
}

uint64_t MonoBasis :: size()
{
    return basis.size();
}

int64_t MonoBasis :: id_of(Tuple &t)
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

std::ostream& operator<< ( std::ostream& os, const MonoBasis& mb )
{
    for( auto it = mb.basis.cbegin(); it != mb.basis.cend(); ++it )
    {
        os << it->id << ": " <<*it << std::endl;
    }
    return os;
}
