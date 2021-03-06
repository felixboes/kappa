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


#ifndef KAPPA_HPP
#define KAPPA_HPP

#include <libhomology/homology.hpp>

#include "cssbasis.hpp"
#include "css.hpp"
#include "factorial.hpp"
#include "generators.hpp"
#include "misc.hpp"
#include "monobasis.hpp"
#include "monocomplex.hpp"
#include "monocochain_field.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.

typedef MonoComplex<ChainComplexQ> MonoComplexQ;
typedef MonoComplex<ChainComplexZm> MonoComplexZm;
typedef MonoComplex<ChainComplexZStorageOnly> MonoComplexZStorageOnly;

typedef MonoCochainField< Q > MonoCochainFieldQ;
typedef MonoCochainField< Zm > MonoCochainFieldZm;

typedef ClusterSpectralSequence<ChainComplexQCSS> ClusterSpectralSequenceQ;
typedef ClusterSpectralSequence<ChainComplexZmCSS> ClusterSpectralSequenceZm;

#endif // KAPPA_HPP
