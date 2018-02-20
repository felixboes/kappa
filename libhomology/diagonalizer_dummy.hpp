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


#ifndef DIAGONALIZER_DUMMY_HPP
#define DIAGONALIZER_DUMMY_HPP

// Description:
//
// This header defines a dummy version of a diagonalizer.
// This is used if we are not at all interested in homology computations, for example when we construct and save differentials.

#include <boost/numeric/ublas/io.hpp>

#include "homology_dummy.hpp"
#include "parallelization.hpp"

/**
 *  This is a dummy class that mimes the functionality of DiagonalozerField, but does nothing.
 */
template <class MatrixType>
class DiagonalizerDummy
{
public:
    DiagonalizerDummy() : transp(false) {}
    void operator() ( MatrixType & ) {}
    uint32_t dfct() {return 0;}
    HomologyDummy::KernT kern() {return HomologyDummy::KernT();}
    uint32_t rank() {return 0;}
    HomologyDummy::TorsT tors() {return HomologyDummy::TorsT();}
    
    void apply_base_changes( MatrixType& differential, const MatrixType& base_changes ){ (void)differential; (void)base_changes; } // Supress 'unused variable' warnings.

    bool transp;
    uint32_t num_working_threads;
    uint32_t num_remaining_threads;
    atomic_uint current_rank;
};

#endif // DIAGONALIZER_DUMMY_HPP
