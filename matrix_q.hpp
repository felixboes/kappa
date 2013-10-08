#ifndef MATRIX_Q_HPP
#define MATRIX_Q_HPP


#include <gmpxx.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

/**
 *  The coefficient ring \f$ \mathbb{Q} \f$.
 */
typedef mpq_class Q;

typedef boost::numeric::ublas::matrix< Q > MatrixQ;
typedef boost::numeric::ublas::identity_matrix< Q > MatrixQIdentity;
typedef boost::numeric::ublas::zero_matrix< Q > MatrixQZero;
typedef boost::numeric::ublas::matrix_slice< MatrixQ > MatrixQSlice;
typedef boost::numeric::ublas::slice slice_q;

#endif // MATRIX_Q_HPP
