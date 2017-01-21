#include "chain_complex.hpp"
#include "chain_complex_impl.ipp"

#include "homology.hpp"

/* Force template instantiation for used types */

template class ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy>;
template class ChainComplex<Q, MatrixQ, DiagonalizerField<MatrixQ>, HomologyField>;
template class ChainComplex<Zm, MatrixZm, DiagonalizerField<MatrixZm>, HomologyField>;
template class ChainComplex<bool, MatrixBool, DiagonalizerBool, HomologyField>;
template class ChainComplex<Q, MatrixQCSS, DiagonalizerField<MatrixQCSS>, HomologyField>;
template class ChainComplex<Zm, MatrixZmCSS, DiagonalizerField<MatrixZmCSS>, HomologyField>;
template class ChainComplex<bool, MatrixBoolCSS, DiagonalizerField<MatrixBoolCSS>, HomologyField>;
