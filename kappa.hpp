#ifndef KAPPA_HPP
#define KAPPA_HPP

#include <homology.hpp>

#include "blockfinder.hpp"
#include "css.hpp"
#include "factorial.hpp"
#include "monocomplex.hpp"
#include "serialization.hpp"
#include "sessionconfig.hpp"
#include "tuple.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.
template class MonoComplex<ChainComplexQ>;
template class MonoComplex<ChainComplexZm>;
template class MonoComplex<ChainComplexZStorageOnly>;

typedef MonoComplex<ChainComplexQ> MonoComplexQ;
typedef MonoComplex<ChainComplexZm> MonoComplexZm;
typedef MonoComplex<ChainComplexZStorageOnly> MonoComplexZStorageOnly;

template class ClusterSpectralSequence<ChainComplexQ>;
template class ClusterSpectralSequence<ChainComplexZm>;
template class ClusterSpectralSequence<ChainComplexZStorageOnly>;

typedef ClusterSpectralSequence<ChainComplexQ> ClusterSpectralSequenceQ;
typedef ClusterSpectralSequence<ChainComplexZm> ClusterSpectralSequenceZm;
typedef ClusterSpectralSequence<ChainComplexZStorageOnly> ClusterSpectralSequenceZStorageOnly;

#endif // KAPPA_HPP
