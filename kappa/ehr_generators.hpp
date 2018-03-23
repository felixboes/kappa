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


#ifndef EHR_GENERATORS_HPP
#define EHR_GENERATORS_HPP

#include <string>

#include "monocomplex.hpp"
#include "monocochain_field.hpp"

enum EhrGenerator
{
    a,
    b,
    c,
    d,
    e,
    f,
    E_1_f,
    E_2_f,
    Qb,
    Qc,
    Qd, ///< This is soley defined for coefficients with 2 = 0.
    Te,
    Eb,
    TEb,
    Q_alpha_inv_c,
    Q_beta_c,
    Q_gamma_c,    
    Q_alpha_inv_d,
    Q_beta_d,
    Q_gamma_d,
    R_a_e,
    R_alpha_inv_beta_c_d,
    R_alpha_inv_gamma_inv_c_d,
    R_alpha_inv_gamma_c_d,
    R_alpha_inv_beta_inv_c_d,
    R_alpha_inv_alpha_inv_c_d,
    R_alpha_inv_beta_c_c,
    R_a_Te, ///< at the moment, it's only defined up to signs. However it is 0 mod 2.
    T_1_f, ///< at the moment, it's only defined up to signs. However it is 0 mod 2.
    T_2_f, ///< at the moment, it's only defined up to signs. However it is 0 mod 2.

    radial_a,
    radial_b,
    radial_Srad_Te,
};

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const EhrGenerator& );

#endif // EHR_GENERATORS_HPP
