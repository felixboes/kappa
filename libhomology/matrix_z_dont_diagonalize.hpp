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


#ifndef MATRIX_Z_DONT_DIAGONALIZE_HPP
#define MATRIX_Z_DONT_DIAGONALIZE_HPP

#include <cinttypes>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

/**
 *  The coefficient ring \f$ \mathbb Z \f$.
 *  Use these matrices only for saving generated differentials.
 *  Do not diagonalize these matrices as the coefficients will probalby overflow.
 *  In order to avoid overflows use matrices with 'gmp' coefficients.
 */

typedef boost::numeric::ublas::matrix< int32_t > MatrixZDontDiagonalize;
typedef boost::numeric::ublas::identity_matrix< int32_t > MatrixZDontDiagonalizeIdentity;
typedef boost::numeric::ublas::zero_matrix< int32_t > MatrixZDontDiagonalizeZero;

#endif // MATRIX_Z_DONT_DIAGONALIZE_HPP
