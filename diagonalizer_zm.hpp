#ifndef DIAGONALIZER_ZM_HPP
#define DIAGONALIZER_ZM_HPP

#include <functional>
#include <iostream>
#include <list>
#include <vector>
#include <stdint.h>

#include "matrix_zm.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *          C_{n-1} \xleftarrow{\partial_{in}} C_n \xleftarrow{\partial_{out}} C_{n-1}
 *  \f]
 *  at position \f$n\f$, it suffices to compute the dimension of the kernel and the image of the respective matrices.
 */

class DiagonalizerZm
{
public:
    typedef std::vector< Zm > TorsionVector;
    DiagonalizerZm(MatrixZm &out_differential, MatrixZm &in_differential );
    uint32_t defect() {return def;}
    TorsionVector torsion_coefficients() {return tors_coefficients;}
    
private:
    
    /**
     *  Performs a row operation to matrix to zeroise the entry (row_2, col) using the entry (row_1, col).
     */
    void row_operation(MatrixZm &matrix, size_t row_1, size_t row_2, size_t col);
    
    /**
     *  Compute the dimension of the image of in.
     *  This done by computing the number of lineary independant columns or rows.
     */
    uint32_t diag_field(MatrixZm& matrix);
    
    MatrixZm &in;
    MatrixZm &out;
    uint32_t def;
    TorsionVector tors_coefficients;
};

#endif // DIAGONALIZER_ZM_HPP
