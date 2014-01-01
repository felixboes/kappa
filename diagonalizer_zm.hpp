#ifndef DIAGONALIZER_ZM_HPP
#define DIAGONALIZER_ZM_HPP

#include <chrono>
#include <cinttypes>
#include <functional>
#include <future>
#include <list>
#include <thread>
#include <vector>

#include <boost/numeric/ublas/io.hpp>

#include "homology_field.hpp"
#include "matrix_zm.hpp"
#include "parallelization.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *          C_{n-1} \xleftarrow{\partial_{out}} C_n \xleftarrow{\partial_{in}} C_{n+1}
 *  \f]
 *  at position \f$n\f$, it suffices to compute the dimension of the kernel and the image of the respective matrices. 
 *  Therefore we need to diagonalize the matrices representing the two differentials (here only for matrices with
 *  coefficients in \f$ \mathbb{F}_p \f$).
 */

class DiagonalizerZm
{
public:
    /**  constructor */
    DiagonalizerZm() {}
    
    /**
    *   Diagonalize a given matrix.
    **/
    void operator() ( MatrixZm &matrix );
    void operator() ( MatrixZm &matrix, atomic_uint & current_rank );
    
    /**
     *  Diagonalize matrix and apply the base change to post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n
     *  \f]
     */
    void operator() ( MatrixZm &post_matrix, MatrixZm &matrix );
    
    /**
     *  Diagonalize matrix and apply the base change to pre_matrix and post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n \xleftarrow{\partial{prematrix}} C_{n+1}
     *  \f]
     */
    void operator() ( MatrixZm &post_matrix, MatrixZm &matrix, MatrixZm &pre_matrix );
    
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
    uint32_t diag_field(MatrixZm& matrix);
    uint32_t diag_field(MatrixZm& matrix, atomic_uint & current_rank);
    
    /**
     *  Performs a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col).
     */
    void row_operation(MatrixZm &matrix, size_t row_1, size_t row_2, size_t col);
    
    MatrixZm in;
    MatrixZm out;
    uint32_t def;
    uint32_t rnk;
};

#endif // DIAGONALIZER_ZM_HPP
