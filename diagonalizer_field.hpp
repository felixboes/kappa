#ifndef DIAGONALIZER_FIELD_HPP
#define DIAGONALIZER_FIELD_HPP

#include <chrono>
#include <cinttypes>
#include <functional>
#include <future>
#include <list>
#include <thread>
#include <vector>

#include "homology_field.hpp"
#include "matrix.hpp"

/**
 *  In order to compute the homology of the chain complex
 *  \f[
 *          C_{n-1} \xleftarrow{\partial_{out}} C_n \xleftarrow{\partial_{in}} C_{n+1}
 *  \f]
 *  at position \f$n\f$, it suffices to compute the dimension of the kernel and the image of the respective matrices. 
 *  Therefore we need to diagonalize the matrices representing the two differentials (here only for matrices with
 *  coefficients in \f$ \mathbb{F}_p \f$).
 */

template < class CoefficientT >
class DiagonalizerField
{
public:
    typedef CoefficientT CoefficientType;
    typedef Matrix<CoefficientType> MatrixType;
    
    /**  constructor */
    DiagonalizerField() {}
    
    /**
    *   Diagonalize a given matrix.
    **/
    void operator() ( MatrixType &matrix );
    void operator() ( MatrixType &matrix, atomic_uint & current_rank );
    
    /**
     *  Diagonalize matrix and apply the base change to post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n
     *  \f]
     */
    void operator() ( MatrixType &post_matrix, MatrixType &matrix );
    
    /**
     *  Diagonalize matrix and apply the base change to pre_matrix and post_matrix.
     *  One should have the following picture in mind.
     *  \f[
     *      C_{n-2} \xleftarrow{\partial_{postmatrix}} C_{n-1} \xleftarrow{\partial_{matrix}} C_n \xleftarrow{\partial{prematrix}} C_{n+1}
     *  \f]
     */
    void operator() ( MatrixType &post_matrix, MatrixType &matrix, MatrixType &pre_matrix );
    
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
    uint32_t diag_field(MatrixType& matrix);
    uint32_t diag_field(MatrixType& matrix, atomic_uint & current_rank);
    
    MatrixType in;
    MatrixType out;
    uint32_t def;
    uint32_t rnk;
};

#include "diagonalizer_field.ipp"

#endif // DIAGONALIZER_FIELD_HPP
