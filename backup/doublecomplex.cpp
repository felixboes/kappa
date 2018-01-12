#include "doublecomplex.hpp"
#include "doublecomplex_impl.ipp"

#define force_template_instantiation( MatrixComplex )\
    template class DoubleComplex< MatrixComplex >;\
    template void doublecomplex_work(DoubleComplex<MatrixComplex> & doublecomplex, DoubleComplexWork & work, const uint32_t p, typename MatrixComplex::MatrixType & differential);

force_template_instantiation(ChainComplexQ)
force_template_instantiation(ChainComplexZm)

#undef force_template_instantiation
