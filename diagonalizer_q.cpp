#include "diagonalizer_q.hpp"

uint32_t DiagonalizerQ::dfct()
{
    return def;
}

HomologyField::KernT DiagonalizerQ::kern()
{
    return def;
}

uint32_t DiagonalizerQ::rank()
{
    return rnk;
}

HomologyField::TorsT DiagonalizerQ::tors()
{
    return rnk;
}

/**
 *  Performs a row operation to matrix to zeroise the entry (row_2, col) using the entry (row_1, col).
 */
void DiagonalizerQ::row_operation(MatrixQ & matrix, size_t row_1, size_t row_2, size_t col)
{
    size_t num_cols = matrix.size2();
    Q c( - ( 1 / matrix(row_1,col) ) * matrix(row_2, col) );
    MatrixQ row_op(2,2);
    
    row_op(0,0) = 1;
    row_op(0,1) = 0;
    row_op(1,0) = c;
    row_op(1,1) = 1;    
   
    MatrixQSlice matrix_slice(matrix, slice_q(row_1, row_2-row_1, 2), slice_q(col,1, num_cols-col) );
    matrix_slice = boost::numeric::ublas::prod( row_op, matrix_slice );
}

/**
 *  Compute the dimension of the image of matrix.
 *  This done by computing the number of lineary independent columns or rows.
 */ 
uint32_t DiagonalizerQ::diag_field(MatrixQ &matrix)
{
    size_t num_rows = matrix.size1();
    size_t num_cols = matrix.size2();
    size_t rank = 0; // Stores the dimension of the subspace spanned by the first 0 \le j < col columns. The dimension is at most of size row.
    
    /// @TODO: Is #rows <= #cols?
    // if( num_rows <= num_cols )
    {
         // We avoid permutations of rows in the Gauss algorithm by performing the following steps:
         //   - Iterate through the columns.
         //   - Iterate through a singly linked list of rows.
         //   - For a given column, remove a row from the list if the common entry of column and row is invertible.
         //   - Apply row operations to all rows beneath the row.
         
        std::list<size_t> rows_to_check;
        
        // Fill the list.
        for( size_t row_1 = 0; row_1 < num_rows; row_1++ )
        {
            rows_to_check.push_back(row_1);
        }
        
        // Iterate through columns. We may stop if the rank is maximal.
        for ( size_t col = 0; col < num_cols && rank < num_rows; col++ )
        {
            // Find first invertible (i.e. non-zero) entry in the remaining rows.
            auto it = rows_to_check.begin();
            for( ; it != rows_to_check.end() && matrix( *it, col ) == 0; ++it )
            {
            }
            
            // Does this row contribute to the rank, i.e. did we find an invertible element?
            if( it != rows_to_check.end() )
            {
                // This row contributes to the rank.
                rank++;
                size_t row_1 = *it;
                // We do not need to check this row in further iterations.
                // After erasing, it will point to the next element in the list.
                it = rows_to_check.erase(it);
                // Use row operations to zeroize the remaining elements of the column.
                for( ; it != rows_to_check.end(); ++it )
                {
                    size_t row_2 = *it;
                    // Assuming that the entry (row_2, col) differs from zero, we perform
                    // a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col). 


                    if( matrix( row_2, col ) != 0 )
                    {
                        row_operation( matrix, row_1, row_2, col );
                    }
                }
            }
        }
    }
    return rank;
}

void DiagonalizerQ::operator() ( MatrixQ &matrix )
{
    rnk = diag_field( matrix );
    def = matrix.size2() - rnk;
}

void DiagonalizerQ::operator() ( MatrixQ &post_matrix, MatrixQ &matrix )
{
    operator ()(matrix);
    std::cerr << "TODO: Additional ops in Diagonalizer::operator() ( MatrixQ &post_matrix, MatrixQ &matrix )" << std::endl;
}

void DiagonalizerQ::operator() ( MatrixQ &post_matrix, MatrixQ &matrix, MatrixQ &pre_matrix )
{
    operator ()(matrix);
    std::cerr << "TODO: Additional ops in Diagonalizer::operator() ( MatrixQ &post_matrix, MatrixQ &matrix, MatrixQ &pre_matrix )" << std::endl;
}

