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


#ifndef HOMOLOGY_DUMMY_HPP
#define HOMOLOGY_DUMMY_HPP

// Description:
//
// This header defines a dummy version of homology.
// This is used if we are not at all interested in homology computations, for example when we construct and save differentials.

#include <iostream>

/**
 *  This is a dummy class that mimes the functionality of HomologyField, but does nothing.
 */
class HomologyDummy
{
public:
    typedef int64_t KernT;
    typedef int64_t TorsT;
    HomologyDummy() {}
    HomologyDummy( int32_t, KernT, TorsT ) {}
    inline void set_kern( int32_t, KernT ) {}
    inline void set_tors( int32_t, TorsT ) {}
    friend std::ostream& operator<< (std::ostream& stream, const HomologyDummy& homol);
};

inline std::ostream& operator<< (std::ostream& stream, const HomologyDummy&)
{
    return stream << "This is a dummy module for homology computations.";
}

#endif // HOMOLOGY_DUMMY_HPP
