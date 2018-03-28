#include "ehr_basis.hpp"
#include "ehr_basis_impl.ipp"


#define force_template_instantiation(TupleT) \
    template class EhrBasis<TupleT>;\
    template std::ostream& operator<< ( std::ostream& stream, const EhrBasis<TupleT> & basis );\
    template EhrBasis<TupleT> load_ehr_basis(const uint32_t g, const uint32_t m, const int32_t p, const bool radial);\

force_template_instantiation(SymGrpTuple)
force_template_instantiation(AltGrpTuple)

#undef force_template_instantiation

template <>
uint EhrBasis<SymGrpTuple> :: add_basis_element_reduced( SymGrpTuple t )
{
    if( t.is_multiple_of_a() )
    {
        return 0;
    }
    else
    {
        return add_basis_element(t);
    }
}