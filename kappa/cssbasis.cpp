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


#include "cssbasis.hpp"

int32_t CSSBasis :: add_basis_element ( Tuple t )
{
    uint32_t num_clusters = t.num_clusters();
    LBasisType& l_basis = basis[num_clusters];
    uint32_t id = l_basis.size();
    t.id = id;
    l_basis.insert( std::move(t) );
    
    return id;
}

int32_t CSSBasis :: size( const int32_t l ) const
{
    return (basis.count(l) != 0 ? basis.at(l).size() : 0 );
}

int32_t CSSBasis :: total_size() const
{
    int32_t size(0);
    for( auto& it : basis )
    {
        size += it.second.size();
    }
    return size;
}

int32_t CSSBasis :: id_of( Tuple& t ) const
{
    for( auto& l_basis_it : basis )
    {
        auto& l_basis = l_basis_it.second;
        auto it = l_basis.find(t);
        if( it != l_basis.end() )
        {
            return it->id;
        }
    }
    return -1;
}

int32_t CSSBasis :: total_id_of( Tuple& t ) const
{
    int64_t basis_offset(0);
    
    for( auto& l_basis_it : basis )
    {
        auto& l_basis = l_basis_it.second;
        auto it = l_basis.find(t);
        if( it != l_basis.end() )
        {
            return basis_offset + it->id;
        }
        basis_offset += this->size(l_basis_it.first);
    }
    return -1;
}

std::ostream& operator<< ( std::ostream& os, const CSSBasis& cb )
{
    for( const auto& it : cb.basis )
    {
        os << "Cluster of size " << it.first << std::endl;
        for( auto& b : it.second )
        {
            os << b.id << ": " << b << std::endl;
        }
    }
    return os;
}
