#ifndef KAPPA_HPP
#define KAPPA_HPP

#include <libhomology/homology.hpp>

#include "blockfinder.hpp"
#include "cssbasis.hpp"
#include "css.hpp"
#include "doublecomplexbasis.hpp"
#include "doublecomplex.hpp"
#include "factorial.hpp"
#include "generators.hpp"
#include "misc.hpp"
#include "monobasis.hpp"
#include "monocomplex.hpp"
#include "monocochain_field.hpp"
#include "sessionconfig.hpp"
#include "operationtester.hpp"
#include "tuple.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.

typedef MonoComplex<ChainComplexQ> MonoComplexQ;
typedef MonoComplex<ChainComplexZm> MonoComplexZm;
typedef MonoComplex<ChainComplexZStorageOnly> MonoComplexZStorageOnly;

typedef MonoCochainField< Q > MonoCochainFieldQ;
typedef MonoCochainField< Zm > MonoCochainFieldZm;

typedef OperationTester< ChainComplexQ,  VectorQ  > OperationTesterQ;
typedef OperationTester< ChainComplexZm, VectorZm > OperationTesterZm;
typedef OperationTesterCSS< ChainComplexQCSS,  VectorQ  > OperationTesterQCSS;
typedef OperationTesterCSS< ChainComplexZmCSS, VectorZm > OperationTesterZmCSS;

typedef ClusterSpectralSequence<ChainComplexQCSS> ClusterSpectralSequenceQ;
typedef ClusterSpectralSequence<ChainComplexZmCSS> ClusterSpectralSequenceZm;

#endif // KAPPA_HPP
