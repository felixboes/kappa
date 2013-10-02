#ifndef DIAGONALIZER_FP_HPP
#define DIAGONALIZER_FP_HPP

#include <list>

#include "matrix_fp.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *          C_{n-1} \xleftarrow{\partial_{post}} C_n \xleftarrow{\partial_{pre}} C_{n-1}
 *  \f]
 *  at position \f$n\f$, it suffices to compute the dimension of the kernel and the image of the respective matrices.
 */

class DiagonalizerFp
{
public:
    DiagonalizerFp( MatrixFp &pre_differential, MatrixFp &post_differential );

private:
    MatrixFp &pre;
    MatrixFp &post;
};

#endif // DIAGONALIZER_FP_HPP
