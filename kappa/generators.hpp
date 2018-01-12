#ifndef GENERATORS_HPP
#define GENERATORS_HPP

#include <string>

#include "monocomplex.hpp"
#include "monocochain_field.hpp"

enum Generator
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
MonoCochainField< CoefficientT > create_cochain( const Generator& );

#endif // GENERATORS_HPP
