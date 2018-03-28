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


#include "ehr_basis.hpp"

template <class TupleT>
EhrBasis<TupleT>::EhrBasis() : basis()
{
}
 
template <class TupleT>
uint32_t EhrBasis<TupleT> :: add_basis_element ( TupleT t )
{
    uint32_t id = basis.size();
    t.id = id;
    basis.insert(std::move(t));
    
    return id;
}

template <class TupleT>
uint EhrBasis<TupleT> :: add_basis_element_reduced( TupleT t )
{
    std::cerr << "Error in 'template <class TupleT> uint EhrBasis<TupleT> :: add_basis_element_reduced( TupleT t )': "
                 "this function should not be called for TupleT unequal to SymGrpTuple! This function was called for "
                 "TupleT t = " << t
              << std::endl;
    return 0;
}

template <class TupleT>
uint64_t EhrBasis<TupleT> :: size() const
{
    return basis.size();
}

template <class TupleT>
int64_t EhrBasis<TupleT> :: id_of(const TupleT &t) const
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

template <class TupleT>
EhrBasis<TupleT> load_ehr_basis(const uint32_t g, const uint32_t m, const int32_t p, const bool radial)
{
    std::string filename =
            "./cache/bases_" + std::string(radial == true ? "radial" : "parallel") + "/" +
            std::to_string(g) + "_" +
            std::to_string(m) + "_" +
            std::to_string(p);

    return load_from_file_bz2< EhrBasis<TupleT> >( filename, false );
}

template <class TupleT>
std::ostream& operator<< ( std::ostream& os, const EhrBasis<TupleT>& mb )
{
    for( const auto it : mb.basis )
    {
        os << it.id << ": " << it << std::endl;
    }
    return os;
}
