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

uint64_t MonoBasis :: size() const
{
    return basis.size();
}

int64_t MonoBasis :: id_of(Tuple &t) const
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

template<>
void update_differential(MatrixField<Zm> & differential,
                         Tuple &           tuple,
                         Tuple &           boundary,
                         int32_t           parity,
                         int8_t            i,
                         int8_t            or_sign,
                         SignConvention &  sign_conv)
{
    int8_t sign_in_differential = sign(parity, i, or_sign, sign_conv);

    if (sign_in_differential == 1)
    {
        differential(tuple.id, boundary.id) += 1;
    }
    else if (sign_in_differential == -1)
    {
        differential(tuple.id, boundary.id) += -1;
    }
}

template<>
void update_differential(MatrixBool &     differential,
                         Tuple &          tuple,
                         Tuple &          boundary,
                         int32_t,
                         int8_t,
                         int8_t,
                         SignConvention &)
{
    differential.add_entry(tuple.id, boundary.id);
}

template<>
void update_differential(MatrixZDontDiagonalize & differential,
                         Tuple &                  tuple,
                         Tuple &                  boundary,
                         int32_t                  parity,
                         int8_t                   i,
                         int8_t                   or_sign,
                         SignConvention &         sign_conv)
{
    int8_t sign_in_differential = sign(parity, i, or_sign, sign_conv);

    if (sign_in_differential == 1)
    {
        differential(tuple.id, boundary.id) += 1;
    }
    else if (sign_in_differential == -1)
    {
        differential(tuple.id, boundary.id) += -1;
    }
}
