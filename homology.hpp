#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

#include "chain_complex.hpp"
#include "clock.hpp"
#include "diagonalizer_dummy.hpp"
#include "diagonalizer_field.hpp"
#include "homology_dummy.hpp"
#include "homology_field.hpp"
#include "matrix.hpp"
#include "matrix_zm.hpp"
#include "matrix_z_dont_diagonalize.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.
//template class ChainComplex<Q, Matrix<Q>, DiagonalizerField<Q>, HomologyField>;
template class ChainComplex<Zm, Matrix<Zm>, DiagonalizerField<Zm>, HomologyField>;
template class ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy>;

//typedef ChainComplex<Q, Matrix<Q>, DiagonalizerField<Q>, HomologyField> ChainComplexQ;
typedef ChainComplex<Zm, Matrix<Zm>, DiagonalizerField<Zm>, HomologyField> ChainComplexZm;
typedef ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy> ChainComplexZStorageOnly;


#endif // HOMOLOGY_HPP
