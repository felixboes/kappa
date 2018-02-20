// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#include "monocomplex.hpp"
#include "monocomplex_impl.ipp"


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

/* Force template instantiation for used types */

// We do not store homchains for Zm right now.
template<>
void MonoComplex< ChainComplexZm > :: homchain(int32_t p, bool homology, int32_t maxdimension)
{
    (void)p;
    (void)homology;
    (void)maxdimension;
}
#define force_template_instantiation(MatrixComplex) \
    template class MonoComplex<MatrixComplex>;\
    template void update_differential(MatrixComplex &differential, const size_t row, const size_t column, const int32_t, const int8_t, const int8_t, const SignConvention &);\
    template void monocomplex_work(MonoComplex<MatrixComplex> & monocomplex, MonocomplexWork & work, const uint32_t p, MatrixComplex::MatrixType & differential);

force_template_instantiation(ChainComplexQ)
force_template_instantiation(ChainComplexZm)
force_template_instantiation(ChainComplexZStorageOnly)

#undef force_template_instantiation\

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

