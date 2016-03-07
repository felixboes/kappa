#ifndef GENERATORS_HPP
#define GENERATORS_HPP

#include <string>

#include "monocomplex.hpp"
#include "operationtester.hpp"

enum Generator
{
    a,
    b,
    c,
    d,
    e,
    Qb,
    Qc,
    Qd, ///< This is soley defined for coefficients with 2 = 0.
    Te,
    Eb,
    TEb,
    Q_alpha_inv_c,
    Q_beta_c,
    Q_gamma_c,    
    Q_alpha_d,
    Q_beta_d,
    Q_gamma_d,
    R_a_e,
    R_alpha_inv_beta_c_d,
    R_alpha_inv_gamma_inv_c_d,
    R_alpha_inv_gamma_c_d,
    R_alpha_inv_beta_inv_c_d,
    R_alpha_inv_alpha_inv_c_d,
    R_alpha_inv_beta_c_c
};

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const Generator& );

#include "generators.ipp"

#endif // GENERATORS_HPP
