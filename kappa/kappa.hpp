#ifndef KAPPA_HPP
#define KAPPA_HPP

#include <libhomology/homology.hpp>

#include "blockfinder.hpp"
#include "css.hpp"
#include "factorial.hpp"
#include "misc.hpp"
#include "monocomplex.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.

typedef MonoComplex<ChainComplexQ> MonoComplexQ;
typedef MonoComplex<ChainComplexZm> MonoComplexZm;
typedef MonoComplex<ChainComplexBool> MonoComplexBool;
typedef MonoComplex<ChainComplexZStorageOnly> MonoComplexZStorageOnly;

#endif // KAPPA_HPP
