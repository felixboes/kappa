#ifndef DIAGONALIZER_FP_HPP
#define DIAGONALIZER_FP_HPP

#include "matrix_fp.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *      \begin{tikzcd}
 *          C_{n-1} \arrow{l}[swap]{\partial_{post}} & C_n \arrow{l}[swap]{\partial_{pre}} & C_{n-1}
 *      \end{tikzcd}
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
