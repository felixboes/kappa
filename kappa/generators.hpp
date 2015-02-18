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
    Qc,
    Qd, ///< This is soley defined for coefficients with 2 = 0.
    Te
};

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const Generator& );

#include "generators.ipp"

#endif // GENERATORS_HPP
