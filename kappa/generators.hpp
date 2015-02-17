#ifndef GENERATORS_HPP
#define GENERATORS_HPP

#include <string>

#include "monocomplex.hpp"

enum Generator
{
    a,
    b
};

template< class CoefficientT >
MonoCochainField< CoefficientT > create_cochain( const Generator& );

#include "generators.ipp"

#endif // GENERATORS_HPP
