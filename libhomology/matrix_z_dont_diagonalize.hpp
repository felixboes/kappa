#ifndef MATRIX_Z_DONT_DIAGONALIZE_HPP
#define MATRIX_Z_DONT_DIAGONALIZE_HPP

#include <cinttypes>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

/**
 *  The coefficient ring \f$ \mathbb{Z} \f$.
 *  Use these matrices only for saving generated differentials.
 *  Do not diagonalize these matrices as the coefficients will probalby overflow.
 *  In order to avoid overflows use matrices with 'gmp' coefficients.
 */

typedef boost::numeric::ublas::matrix< int32_t > MatrixZDontDiagonalize;
typedef boost::numeric::ublas::identity_matrix< int32_t > MatrixZDontDiagonalizeIdentity;
typedef boost::numeric::ublas::zero_matrix< int32_t > MatrixZDontDiagonalizeZero;

#endif // MATRIX_Z_DONT_DIAGONALIZE_HPP
