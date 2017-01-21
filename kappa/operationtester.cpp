#include "operationtester.hpp"
#include "operationtester_impl.ipp"

#define force_template_instantiation( CoefficientType, MatrixComplex, VectorType  )\
    template class OperationTester< MatrixComplex,  VectorType  >;\
    template VectorType kappa_dual( const CoefficientType& c, const Tuple& t, const MonoBasis& b );\
    template void compute_and_add_kappa_dual_rec( const CoefficientType& c, const Tuple& t, const MonoBasis& b, VectorType& v, const std::vector<size_t> s, const size_t i );


force_template_instantiation(Q, ChainComplexQ, VectorQ)
force_template_instantiation(Zm, ChainComplexZm, VectorZm)

#undef force_template_instantiation

typedef OperationTesterCSS< ChainComplexQCSS,  VectorQ  > OperationTesterQCSS;
typedef OperationTesterCSS< ChainComplexZmCSS, VectorZm > OperationTesterZmCSS;
