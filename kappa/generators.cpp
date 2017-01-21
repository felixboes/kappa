#include "generators.hpp"
#include "generators_impl.ipp"

#include "kappa.hpp"

/* Force template instantiation for used types */

template MonoCochainField< Q > create_cochain( const Generator& );
template MonoCochainField< Zm > create_cochain( const Generator& );
