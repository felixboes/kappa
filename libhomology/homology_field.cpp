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


#include "homology_field.hpp"

HomologyField::HomologyField()
{
}

HomologyField::HomologyField( int32_t n, KernT ker, TorsT tor )
{
    kern[n] = ker;
    tors[n] = tor;
}

void HomologyField::set_kern( int32_t n, KernT k )
{
    kern[n] = k;
}

void HomologyField::set_tors( int32_t n, TorsT t )
{
    tors[n] = t;
}

HomologyField::KernT HomologyField::get_kern( int32_t n ) const
{
    if( kern.count(n) == 0 )
    {
        return KernT();;
    }
    else
    {
        return kern.at(n);
    }
}

HomologyField::TorsT HomologyField::get_tors( int32_t n ) const
{
    if( tors.count(n) == 0 )
    {
        return TorsT();
    }
    else
    {
        return tors.at(n);
    }
}

void HomologyField::erase_kern( int32_t n )
{
    kern.erase(n);
}

void HomologyField::erase_tors( int32_t n )
{
    tors.erase(n);
}

std::ostream &operator << ( std::ostream& os, const HomologyField& homol )
{
    for( const auto& it : homol.kern )
    {
        auto & p = it.first;
        auto & kern = it.second;
        
        os << "Homology module H_" << p << std::endl;
        os << std::setfill('-') << std::setw(35) << "-" << std::setfill(' ') << std::endl;
        os << "Dimension = " << kern - homol.get_tors(p) << std::endl;
        os << std::endl;
    }
    return os;
}
