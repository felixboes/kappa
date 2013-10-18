#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

#include "chain_complex.hpp"
#include "diagonalizer_q.hpp"
#include "diagonalizer_zm.hpp"
#include "homology_field.hpp"
#include "matrix_q.hpp"
#include "matrix_zm.hpp"

// In order to use chain complexes with rational and Zm coefficients in other projects,
// we have to use instanciate the templates explicitly.
template class ChainComplex<Q, MatrixQ, DiagonalizerQ, HomologyField>;
template class ChainComplex<Zm, MatrixZm, DiagonalizerZm, HomologyField>;

typedef ChainComplex<Q, MatrixQ, DiagonalizerQ, HomologyField> ChainComplexQ;
typedef ChainComplex<Zm, MatrixZm, DiagonalizerZm, HomologyField> ChainComplexZm;


#endif // HOMOLOGY_HPP
