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
    TEb
};

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const Generator& );

#include "generators.ipp"

#endif // GENERATORS_HPP
