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


#ifndef HOMOLOGY_FIELD_HPP
#define HOMOLOGY_FIELD_HPP

// Description:
//
// This header defines the homology for chaincomplexes over a field.
// Till now we are only interested in the dimensions of kernels and images.

#include <cassert>
#include <cinttypes>
#include <iomanip>
#include <iostream>
#include <map>

/**
 * Class implementing the homology of a chain complex with coefficients in a field.
 * Over a field, a homology module is characterized by the defect of the
 * outgoing differential and the dimension of the rank of the ingoing differential.
 */
class HomologyField
{
public:
    /// type of the kernel
    typedef int32_t KernT;

    /// type of the torsion
    typedef int32_t TorsT;
    
    /// trivial constructor
    HomologyField();

    /// constructor initializing a single homology module
    HomologyField( int32_t, KernT, TorsT );
    
    /// sets the kernel
    void set_kern( int32_t n, KernT k );

    /// sets the torsion
    void set_tors( int32_t n, TorsT t );
    
    /// gets the kernel
    KernT get_kern( int32_t ) const;
    
    /// sets the torsion
    TorsT get_tors( int32_t ) const;
    
    /// erase the kernel
    void erase_kern( int32_t n );
    
    /// erase the torsion
    void erase_tors( int32_t n );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const HomologyField& homol);

protected:
    /// stores the kernel for each possibly non-zero homology module
    std::map< int32_t, int32_t > kern;

    /// stores the torsion for each possibly non-zero homology module
    std::map< int32_t, int32_t > tors;
};

std::ostream& operator<< (std::ostream& stream, const HomologyField& homol);

#endif // HOMOLOGY_FIELD_HPP
