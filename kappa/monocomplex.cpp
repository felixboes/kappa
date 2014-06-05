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
void update_differential(MatrixBool &           differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          ,
                         const int8_t           ,
                         const int8_t           ,
                         const SignConvention & )
{
    differential.add_entry(row, column);
}

template<>
void update_differential(MatrixBoolCSS &        differential,
                         const size_t           row,
                         const size_t           column,
                         const int32_t          ,
                         const int8_t           ,
                         const int8_t           ,
                         const SignConvention & )
{
    differential.add_entry(row, column);
}

int32_t sign(const int32_t          parity,
             const int8_t           i,
             const int8_t           or_sign,
             const SignConvention & sign_conv )
{
    if ( sign_conv == no_signs)
    {
        return 1;
    }
    if( sign_conv == all_signs )
    {
        int32_t actual_parity = (parity + i) % 2;
        if ( or_sign == -1 )
        {
            actual_parity = (actual_parity + 1) % 2;
        }
        //std::cout << it << " " << i << ": The d^hor_i boundary of " << current_basis << ". This is " << boundary << std::endl;
        //std::cout << it.id << "->" << boundary.id << " in " << "M_{" << basis_complex[p-1].size() << "," << basis_complex[p].size() << "} parity=" << actual_parity << std::endl;
        //std::cout << std::endl;
        if ( actual_parity == 0 )
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
    else if( sign_conv == no_orientation_sign )
    {
        if ( (parity + i) % 2 == 0 )
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    return 0;
}

