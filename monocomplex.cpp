#include "monocomplex.hpp"

/*
 *
 *   MonoBasis
 *
 */
void MonoBasis :: add_basis_element (Tuple& t)
{
    basis.push_back(t);
}


uint64_t MonoBasis :: size()
{
    return basis.size();
}

std::ostream& operator<< ( std::ostream& os, const MonoBasis& mb )
{
    for( auto it = mb.basis.cbegin(); it != mb.basis.cend(); ++it )
    {
        os << *it << std::endl;
    }
    return os;
}
