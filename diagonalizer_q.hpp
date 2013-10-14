#ifndef DIAGONALIZER_Q_HPP
#define DIAGONALIZER_Q_HPP

#include <cinttypes>
#include <iostream>
#include <list>
#include <vector>

#include "homology_field.hpp"
#include "matrix_q.hpp"

class DiagonalizerQ
{
public:
    /**  constructor */
    DiagonalizerQ() {}
    
    /**
    *   Diagonalize a given matrix.
    **/
    void operator() ( MatrixQ &matrix );
    
    /**
     *  Diagonalize matrix and apply the base change to post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n
     *  \f]
     */
    void operator() ( MatrixQ &post_matrix, MatrixQ &matrix );
    
    /**
     *  Diagonalize matrix and apply the base change to pre_matrix and post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n \xleftarrow{\partial{prematrix}} C_{n+1}
     *  \f]
     */
    void operator() ( MatrixQ &post_matrix, MatrixQ &matrix, MatrixQ &pre_matrix );
    
    /**  @return defect of the matrix */
    uint32_t dfct();   
    /**  @return defect of the matrix */
    HomologyField::KernT kern();
    /**  @return rank of the matrix */
    uint32_t rank();
    /**  @return rank of the matrix */
    HomologyField::TorsT tors();
    
private:
    /**
     *  @return rank of matrix
     *  The matrix is diagonalized via Gauss to compute the number of linearly independant columns or rows.
     */
    uint32_t diag_field(MatrixQ& matrix);
    
    /**
     *  Performs a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col).
     */
    void row_operation(MatrixQ &matrix, size_t row_1, size_t row_2, size_t col);
    
    MatrixQ in;
    MatrixQ out;
    uint32_t def;
    uint32_t rnk;
};


#endif // DIAGONALIZER_Q_HPP
