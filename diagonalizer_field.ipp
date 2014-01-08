#include "diagonalizer_field.hpp"

template< class CoefficientT >
uint32_t DiagonalizerField< CoefficientT >::dfct()
{
    return def;
}

template< class CoefficientT >
HomologyField::KernT DiagonalizerField< CoefficientT >::kern()
{
    return def;
}

template< class CoefficientT >
uint32_t DiagonalizerField< CoefficientT >::rank()
{
    return rnk;
}

template< class CoefficientT >
HomologyField::TorsT DiagonalizerField< CoefficientT >::tors()
{
    return rnk;
}

/**
 *  Compute the dimension of the image of matrix.
 *  This done by computing the number of lineary independent columns or rows.
 */ 
template< class CoefficientT >
uint32_t DiagonalizerField< CoefficientT >::diag_field(MatrixType &matrix, atomic_uint& current_rank)
{
    size_t num_rows = matrix.size1();
    size_t num_cols = matrix.size2();
    current_rank = 0; // Stores the dimension of the subspace spanned by the first 0 \le j < col columns. The dimension is at most of size row.
    
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
        for( size_t row_1 = 0; row_1 < num_rows; ++row_1 )
        {
            rows_to_check.push_back(row_1);
        }
        
        size_t col = 0;
        size_t row_1 = 0;
        
        // Iterate through columns. We may stop if the rank is maximal.
        for (; col < num_cols && current_rank < num_rows; ++col )
        {
            // Find first invertible (i.e. non-zero) entry in the remaining rows.
            auto it = rows_to_check.begin();
            for( ; it != rows_to_check.end() && matrix( *it, col ) == CoefficientType(0); ++it )
            {
            }
            
            // Does this row contribute to the rank, i.e. did we find an invertible element?
            if( it != rows_to_check.end() )
            {
                // This row contributes to the rank.
                current_rank++;
                row_1 = *it;
                // We do not need to check this row in further iterations.
                // After erasing, it will point to the next element in the list.
                it = rows_to_check.erase(it);
                // Use row operations to zeroize the remaining elements of the column.
                for( size_t relative_position=0 ; it != rows_to_check.end(); ++it, ++relative_position )
                {
                    size_t row_2 = *it;
                    // Assuming that the entry (row_2, col) differs from zero, we perform
                    // a row operation to matrix to zeroize the entry (row_2, col) using the entry (row_1, col). 
                    if( matrix( row_2, col ) != CoefficientType(0) )
                    {
                        matrix.row_operation( row_1, row_2, col );
                    }
                }
            }
        }
    }
    return current_rank;
}

template< class CoefficientT >
void DiagonalizerField< CoefficientT >::operator() ( MatrixType &matrix )
{
    atomic_uint current_rank;
    this->operator ()(matrix, current_rank);
}

template< class CoefficientT >
void DiagonalizerField< CoefficientT >::operator() ( MatrixType &matrix, atomic_uint & current_rank )
{
    rnk = diag_field( matrix, current_rank );
    def = matrix.size2() - rnk;
}

template< class CoefficientT >
void DiagonalizerField< CoefficientT >::operator() ( MatrixType &post_matrix, MatrixType &matrix )
{
    operator ()(matrix);
    std::cerr << "TODO: Additional ops in DiagonalizerField::operator() ( MatrixType &post_matrix, MatrixType &matrix )" << std::endl;
}

template< class CoefficientT >
void DiagonalizerField< CoefficientT >::operator() ( MatrixType &post_matrix, MatrixType &matrix, MatrixType &pre_matrix )
{
    operator ()(matrix);
    std::cerr << "TODO: Additional ops in DiagonalizerField::operator() ( MatrixType &post_matrix, MatrixType &matrix, MatrixType &pre_matrix )" << std::endl;
}
