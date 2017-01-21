#include "monocochain_field.hpp"
#include "monocochain_field_impl.ipp"

#define force_template_instantiation(Coeff) \
    template class MonoCochainField<Coeff>;\
    template std::ostream& operator<< ( std::ostream& stream, const MonoCochainField<Coeff> & cochain );\
    template MonoCochainField< Coeff > operator*( const MonoCochainField< Coeff >&  x, const MonoCochainField< Coeff >& y );

force_template_instantiation(Q)
force_template_instantiation(Zm)

#undef force_template_instantiation
