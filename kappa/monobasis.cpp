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


#include "monobasis.hpp"

MonoBasis::MonoBasis() : basis()
{
}
 
uint32_t MonoBasis :: add_basis_element ( SymGrpTuple t )
{
    uint32_t id = basis.size();
    t.id = id;
    basis.insert(std::move(t));
    
    return id;
}

uint MonoBasis :: add_basis_element_reduced( SymGrpTuple t )
{
    if( t.is_multiple_of_a() )
    {
        return 0;
    }
    else
    {
        return add_basis_element(t);
    }
}

uint64_t MonoBasis :: size() const
{
    return basis.size();
}

int64_t MonoBasis :: id_of(const SymGrpTuple &t) const
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

MonoBasis load_mono_basis( const uint32_t g, const uint32_t m, const int32_t p, const bool radial )
{
    std::string filename =
            "./cache/bases_" + std::string(radial == true ? "radial" : "parallel") + "/" +
            std::to_string(g) + "_" +
            std::to_string(m) + "_" +
            std::to_string(p);

    return load_from_file_bz2< MonoBasis >( filename, false );
}

std::ostream& operator<< ( std::ostream& os, const MonoBasis& mb )
{
    for( const auto it : mb.basis )
    {
        os << it.id << ": " << it << std::endl;
    }
    return os;
}
