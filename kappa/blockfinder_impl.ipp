#include "blockfinder.hpp"


/**
 * Determines all columns c with matrix(row, c) non-zero that have not been visited before,
 * puts them into stack_cols and marks them as visited.
 */
template < class MatrixT >
void add_incident_cols(int                 row,
                       MatrixT &           matrix,
                       std::stack< int > & stack_cols,
                       std::vector<bool> & visited_cols)
{
    int num_cols = matrix.size2();
    for ( int col = 0; col < num_cols; ++col )
    {
        if ( (bool)(matrix(row, col) ) && (visited_cols[col] == false) )
        {
            stack_cols.push(col);
            visited_cols[col] = true;
        }
    }
}

/**
 * Determines all rows r with matrix(r, col) non-zero that have not been visited before,
 * puts them into stack_rows and marks them as visited.
 */
template< class MatrixT >
void add_incident_rows(int                 col,
                       MatrixT &           matrix,
                       std::stack< int > & stack_rows,
                       std::vector<bool> & visited_rows)
{
    int num_rows = matrix.size1();
    for ( int row = 0; row < num_rows; ++row )
    {
        if ( (bool)(matrix(row, col) ) && (visited_rows[row] == false) )
        {
            stack_rows.push(row);
            visited_rows[row] = true;
        }
    }
}

/**
 * Determines the block of the matrix to which the row initial_row belongs.
 * Thereby we store in the arrays visited_rows and visited_cols whether we 
 * have already considered a row or a column.
 */
template< class MatrixT >
static Block find_block(int initial_row,
                        MatrixT & matrix,
                        std::vector<bool> & visited_rows,
                        std::vector<bool> & visited_cols)
{
    Component comp_rows;
    Component comp_cols;
    
    // At the end of this function, the block of initial_row will consist of the rows comp_rows
    // and the columns comp_cols.
    // We start at the row initial_row and put all columns c for which the entry (initial_row, c)
    // is non-zero into the stack_cols. Whenever there is a row r_1 left in stack_rows (resp. a column c_1
    // in stack_columns), we put all columns c_2 with matrix(r_1, c_2) non-zero (resp. all rows r_2 with 
    // matrix(r_2, c_1) non-zero) onto the corresponding stack and continue. Thereby we ignore rows or
    // columns which we have visited before since they were already considered in this block or another.
    
    std::stack< int > stack_rows;
    std::stack< int > stack_cols;

    stack_rows.push(initial_row);

    while ( stack_rows.size() != 0 or stack_cols.size() != 0 )
    {
        if ( stack_rows.size() != 0 )
        {
            int row = stack_rows.top();
            add_incident_cols(row, matrix, stack_cols, visited_cols);
            comp_rows.push_back(row);
            stack_rows.pop();
        }
        else
        {
            int col = stack_cols.top();
            add_incident_rows(col, matrix, stack_rows, visited_rows);
            comp_cols.push_back(col);
            stack_cols.pop();
        }
    }
    return std::make_pair(comp_rows, comp_cols);
}

template< class MatrixT >
BlockFinder<MatrixT>::BlockFinder(MatrixType & matrix)
{
    int num_rows = matrix.size1();
    int num_cols = matrix.size2();
    std::vector<bool> visited_rows(num_rows, false); ///< stores whether a row has been visited
    std::vector<bool> visited_cols(num_cols, false); ///< stores whether a column has been visited
    // We iterate through the rows of the matrix. For a given row r, we check whether it has been visited
    // before, and if not, we determine the block to which r belongs.
    for ( int row = 0; row < num_rows; ++row )
    {
        if ( visited_rows[row] == false )
        {
            _block_part.push_back(find_block(row, matrix, visited_rows, visited_cols));
        }
    }   
}

template< class MatrixT >
BlockPartition & BlockFinder<MatrixT>::block_partition()
{
    return _block_part;
}

template< class MatrixT >
int BlockFinder<MatrixT>::num_blocks() const
{
    return _block_part.size();
}

template< class MatrixT >
int BlockFinder<MatrixT>::num_non_zero_blocks() const
{
    int num = 0;
    for (auto & block : _block_part)
    {
        // A block is a zero-row iff it contains a single row and no column. 
        if (not (block.first.size() == 1 && block.second.size() == 0))
        {
            ++num;
        }
    }
    return num;
}
