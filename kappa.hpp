#ifndef KAPPA_HPP
#define KAPPA_HPP

#include <homology.hpp>

#include "factorial.hpp"
#include "monocomplex.hpp"
#include "tuple.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.
template class MonoComplex<ChainComplexQ>;
template class MonoComplex<ChainComplexZm>;

typedef MonoComplex<ChainComplexQ> MonoComplexQ;
typedef MonoComplex<ChainComplexZm> MonoComplexZm;

#endif // KAPPA_HPP
