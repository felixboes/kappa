#ifndef DIAGONALIZER_ZM_HPP
#define DIAGONALIZER_ZM_HPP

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
 *  Therefore we need to diagonalize the matrices representing the two differentials (here only for matrices with
 *  coefficients in Z/m).
 */

class DiagonalizerZm
{
public:
    typedef std::vector< Zm > TorsionVector;  
    
    /**  constructor */
    DiagonalizerZm(MatrixZm &out_differential, MatrixZm &in_differential );
    
    /**  \return defect of the matrix out */
    uint32_t defect() {return def;}
    
    /**  \return vector of torsion coefficients in the homology of C_n (with the above notation) */
    TorsionVector torsion_coefficients() {return tors_coefficients;}
    
private:
    
    /**
     *  Performs a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col).
     */
    void row_operation(MatrixZm &matrix, size_t row_1, size_t row_2, size_t col);
    
    /**
     *  \return rank of matrix
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     */
    uint32_t diag_field(MatrixZm& matrix);
    
    MatrixZm &in;
    MatrixZm &out;
    uint32_t def;
    TorsionVector tors_coefficients;
};

#endif // DIAGONALIZER_ZM_HPP
