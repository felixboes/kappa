// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

// Description:
//
// This header should be used as in include in other projects.
// It defines chaincomplexes with coefficients in Q and Zm and offers homology computations.
// It defines chaincomplexes with coefficients in Z but denies homology computations.

#include "chain_complex.hpp"
#include "clock.hpp"
#include "condition.hpp"
#include "field_coefficients.hpp"
#include "diagonalizer_dummy.hpp"
#include "diagonalizer_field.hpp"
#include "homology_dummy.hpp"
#include "homology_field.hpp"
#include "matrix_field.hpp"
#include "matrix_z_dont_diagonalize.hpp"
#include "parallelization.hpp"
#include "serialization.hpp"
#include "thread.hpp"
#include "vector_field.hpp"


// Coefficients
typedef ZmBase<> Zm;

// Vectors
typedef VectorField<Q> VectorQ;     ///< This defines Vectors with \f$\mathbb Q\f$ coefficients.
typedef VectorField<Zm> VectorZm;   ///< This defines Vectors with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

// Matrices
typedef MatrixField<Q> MatrixQ;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixField< Zm >MatrixZm;  ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.
typedef MatrixFieldCSS<Q> MatrixQCSS;     ///< This defines Matrices with \f$\mathbb Q\f$ coefficients.
typedef MatrixFieldCSS< Zm > MatrixZmCSS; ///< This defines Matrices with \f$\mathbb Z/ m\mathbb Zf$ coefficients.

// Diagonalizers
typedef DiagonalizerField<MatrixQ> DiagonalizerQ;             ///< This defines Diagonalizer with \f$\mathbb Q\f$ coefficients.
typedef DiagonalizerField<MatrixZm> DiagonalizerZm;           ///< This defines Diagonalizer with \f$\mathbb Z_m\f$ coefficients.
typedef DiagonalizerField<MatrixBool> DiagonalizerBool;       ///< This defines Diagonalizer with \f$\Z_2\f$ coefficients.
typedef DiagonalizerField<MatrixQCSS> DiagonalizerQCSS;       ///< This defines Diagonalizer with \f$\mathbb Q\f$ coefficients and CSS matrices.
typedef DiagonalizerField<MatrixZmCSS> DiagonalizerZmCSS;     ///< This defines Diagonalizer with \f$\mathbb Z_m\f$ coefficients and CSS matrices.
typedef DiagonalizerField<MatrixBoolCSS> DiagonalizerBoolCSS; ///< This defines Diagonalizer with \f$\Z_2\f$ coefficients and CSS matrices.

//Chaincomplexes
typedef ChainComplex<int32_t, MatrixZDontDiagonalize, DiagonalizerDummy<MatrixZDontDiagonalize>, HomologyDummy> ChainComplexZStorageOnly;
typedef ChainComplex<Q, MatrixQ, DiagonalizerField<MatrixQ>, HomologyField> ChainComplexQ;
typedef ChainComplex<Zm, MatrixZm, DiagonalizerField<MatrixZm>, HomologyField> ChainComplexZm;
typedef ChainComplex<bool, MatrixBool, DiagonalizerBool, HomologyField> ChainComplexBool;
typedef ChainComplex<Q, MatrixQCSS, DiagonalizerField<MatrixQCSS>, HomologyField> ChainComplexQCSS;
typedef ChainComplex<Zm, MatrixZmCSS, DiagonalizerField<MatrixZmCSS>, HomologyField> ChainComplexZmCSS;
typedef ChainComplex<bool, MatrixBoolCSS, DiagonalizerField<MatrixBoolCSS>, HomologyField> ChainComplexBoolCSS;

#endif // HOMOLOGY_HPP
