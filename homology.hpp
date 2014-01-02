#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

#include "chain_complex.hpp"
#include "clock.hpp"
#include "diagonalizer_dummy.hpp"
#include "diagonalizer_q.hpp"
#include "diagonalizer_zm.hpp"
#include "homology_dummy.hpp"
#include "homology_field.hpp"
#include "matrix_z_dont_diagonalize.hpp"
#include "matrix_q.hpp"
#include "matrix_zm.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.
template class ChainComplex<Q, MatrixQ, DiagonalizerQ, HomologyField>;
template class ChainComplex<Zm, MatrixZm, DiagonalizerZm, HomologyField>;
template class ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy>;

typedef ChainComplex<Q, MatrixQ, DiagonalizerQ, HomologyField> ChainComplexQ;
typedef ChainComplex<Zm, MatrixZm, DiagonalizerZm, HomologyField> ChainComplexZm;
typedef ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy> ChainComplexZStorageOnly;


#endif // HOMOLOGY_HPP
