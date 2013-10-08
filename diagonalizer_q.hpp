#ifndef DIAGONALIZER_Q_HPP
#define DIAGONALIZER_Q_HPP

#include <iostream>
#include <list>
#include <vector>
#include <stdint.h>

#include "matrix_q.hpp"

class DiagonalizerQ
{
public:
    /**  constructor */
    DiagonalizerQ() {}
    
    /**  constructor */
    DiagonalizerQ(MatrixQ &out_differential, MatrixQ &in_differential );
    uint32_t defect () { return def; }
    uint32_t torsion() { return tor; }
    /**
     *  Compute the dimension of the image of matrix.
     *  This done by computing the number of lineary independant columns or rows.
     */
    uint32_t diag_field(MatrixQ& matrix);

private:
    
    /**
     *  Performs a row operation to matrix to zeroise the entry (row_2, col) using the entry (row_1, col).
     */
    void row_operation(MatrixQ &matrix, size_t row_1, size_t row_2, size_t col);
    
    uint32_t def;
    uint32_t tor;
    MatrixQ in;
    MatrixQ out;
};


#endif // DIAGONALIZER_Q_HPP
