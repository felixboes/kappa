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


#include "diagonalizer_field.hpp"
#include "diagonalizer_field_impl.ipp"

#include "homology.hpp"

// For bool versions of the matrices, we do not have an efficient base change algorithm yet.
template<>
void DiagonalizerField< MatrixBool >::apply_base_changes( MatrixType& differential, const MatrixType& base_changes )
{
    (void)differential;
    (void)base_changes;
}

template<>
void DiagonalizerField< MatrixBoolCSS >::apply_base_changes( MatrixType& differential, const MatrixType& base_changes )
{
    (void)differential;
    (void)base_changes;
}


#define force_template_instantiation(MatrixType)\
    template class DiagonalizerField<MatrixType>;

force_template_instantiation(MatrixQ)
force_template_instantiation(MatrixZm)
force_template_instantiation(MatrixBool)
force_template_instantiation(MatrixQCSS)
force_template_instantiation(MatrixZmCSS)
force_template_instantiation(MatrixBoolCSS)

#undef force_template_instantiation

