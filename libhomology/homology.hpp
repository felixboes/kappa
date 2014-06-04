#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

// Description:
//
// This header should be used as in include in other projects.
// It defines chaincomplexes with coefficients in Q and Zm and offers homology computations.
// It defines chaincomplexes with coefficients in Z but denies homology computations.

#include "chain_complex.hpp"
#include "clock.hpp"
#include "field_coefficients.hpp"
#include "diagonalizer_dummy.hpp"
#include "diagonalizer_field.hpp"
#include "homology_dummy.hpp"
#include "homology_field.hpp"
#include "matrix_field.hpp"
#include "matrix_z_dont_diagonalize.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.
template class ChainComplex<Q, MatrixQ, DiagonalizerField<MatrixQ>, HomologyField>;
template class ChainComplex<Zm, MatrixZm, DiagonalizerField<MatrixZm>, HomologyField>;
template class ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy>;

typedef ChainComplex<Q, MatrixQ, DiagonalizerField<MatrixQ>, HomologyField> ChainComplexQ;
typedef ChainComplex<Zm, MatrixZm, DiagonalizerField<MatrixZm>, HomologyField> ChainComplexZm;
typedef ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy> ChainComplexZStorageOnly;
typedef ChainComplex<bool, MatrixBool, DiagonalizerField<MatrixBool>, HomologyField> ChainComplexBool;
typedef ChainComplex<Q, MatrixCSSQ, DiagonalizerField<MatrixCSSQ>, HomologyField> ChainComplexQCSS;
typedef ChainComplex<Zm, MatrixCSSZm, DiagonalizerField<MatrixCSSZm>, HomologyField> ChainComplexZmCSS;


#endif // HOMOLOGY_HPP
