#ifndef DIAGONALIZER_ZM_HPP
#define DIAGONALIZER_ZM_HPP

#include <cinttypes>
#include <functional>
#include <iostream>
#include <list>
#include <vector>
#include <boost/numeric/ublas/io.hpp>

#include "matrix_zm.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *          C_{n-1} \xleftarrow{\partial_{in}} C_n \xleftarrow{\partial_{out}} C_{n-1}
 *  \f]
 *  at position \f$n\f$, it suffices to compute the dimension of the kernel and the image of the respective matrices. 
 *  Therefore we need to diagonalize the matrices representing the two differentials (here only for matrices with
 *  coefficients in Z/m).
 */

class DiagonalizerZm
{
public:
    /**  constructor */
    DiagonalizerZm() {}
    
    /**  constructor */
    DiagonalizerZm(MatrixZm &out_differential, MatrixZm &in_differential );
    
    /**  @return defect of the matrix out */
    uint32_t defect();   
    /**  @return defect of the matrix out */
    uint32_t kern();
    /**  @return rank of the matrix in */
    uint32_t torsion();
    
    /**
     *  @return rank of matrix
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     */
    uint32_t diag_field(MatrixZm& matrix);
    
private:
    /**
     *  Performs a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col).
     */
    void row_operation(MatrixZm &matrix, size_t row_1, size_t row_2, size_t col);
    
    MatrixZm in;
    MatrixZm out;
    uint32_t def;
    uint32_t tor;
};

#endif // DIAGONALIZER_ZM_HPP
