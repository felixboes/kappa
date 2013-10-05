#include "diagonalizer_q.hpp"

DiagonalizerQ::DiagonalizerQ(MatrixQ &out_differential, MatrixQ &in_differential ) :
    out(out_differential),
    in(in_differential),
    def(0),
    tor(0)
{
    // Defect equals #cols - rank
    def = out.size2() - diag_field(out);
    tor = diag_field(in);
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
 *  This done by computing the number of lineary independant columns or rows.
 */ 
uint32_t DiagonalizerQ::diag_field(MatrixQ &matrix)
{
    size_t num_rows = matrix.size1();
    size_t num_cols = matrix.size2();
    size_t rank = 0; // Stores the dimension of the subspace spanned by the first 0 \le j < col_c columns. The dimension is at most of size row_c.
    
    /// @TODO: Is #rows <= #cols?
    // if( num_rows <= num_cols )
    {
        // Wir gehen davon aus, dass die Anzahl der Zeilen echt kleiner als die Anzahl der Spalten ist.
        //
        // Um Dimension des Kerns und des Bildes zu berechnen, reicht es aus,
        // die Anzahl der linear unabhängigen Zeilen zu bestimmen.
        //
        // Wir können im Stufenform-Algorithmus Vertauschungen zu vermeiden:
        // - Iteriere durch die Spalten.
        //   - Iteriere durch einfach verkettete Liste der Zeilen.
        //   - Lösche Eintrag, falls es in dieser Zeile (bei gegebener Spalte) ein invertierbares Element gibt.
        //   - Mache Zeilenop auf alle vorhandenen Einträge.
        
        std::list<size_t> rows_to_check;
        
        // Fill the list.
        for( size_t row_1 = 0; row_1 < num_rows; row_1++ )
        {
            rows_to_check.push_back(row_1);
        }
        
        // Iterate through columns. We may stop if the rank is maximal.
        for ( size_t col = 0; col < num_cols && rank < num_rows; col++ )
        {
            // Find first non-zero entry in the remaining rows.
            auto it = rows_to_check.begin();
            for( ; it != rows_to_check.end() && matrix( *it, col ) == 0; ++it )
            {
            }
            
            // Does this row contribute to the rank, i.e. did we find an invertable element?
            if( it != rows_to_check.end() )
            {
                // This row contributes to the rank.
                rank++;
                size_t row_1 = *it;
                // We do not need to check this row in further iterations.
                // After erasing, it will point to the next element in the list.
                it = rows_to_check.erase(it);
                // Use row operations to zeroise the remaining elements of the column.
                for( ; it != rows_to_check.end(); ++it )
                {
                    size_t row_2 = *it;
                    // Assume that the entry (row_2, col) differs from zero.
                    // Performs a row operation to matrix to zeroise the entry (row_2, col) 
                    // using the entry (row_1, col). 
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
