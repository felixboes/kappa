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


#include "ehr_cochain_field.hpp"
#include "ehr_cochain_field_impl.ipp"

#define force_template_instantiation(Coeff) \
    template class EhrCochainField<Coeff>;\
    template std::ostream& operator<< ( std::ostream& stream, const EhrCochainField<Coeff> & cochain );\
    template EhrCochainField< Coeff > operator*( const EhrCochainField< Coeff >&  x, const EhrCochainField< Coeff >& y );

force_template_instantiation(Q)
force_template_instantiation(Zm)

#undef force_template_instantiation
