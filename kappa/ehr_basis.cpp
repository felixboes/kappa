#include "ehr_basis.hpp"
#include "ehr_basis_impl.ipp"


#define force_template_instantiation(TupleT) \
    template class EhrBasis<TupleT>;\
    template std::ostream& operator<< ( std::ostream& stream, const EhrBasis<TupleT> & basis );\
    template EhrBasis<TupleT> load_ehr_basis(const uint32_t g, const uint32_t m, const int32_t p, const bool radial);\

force_template_instantiation(SymGrpTuple)

#undef force_template_instantiation