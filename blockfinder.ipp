#include "blockfinder.hpp"

BlockFinder::BlockFinder()
{
}

void add_incident_cols(int                 row,
                       MatrixT &           matrix,
                       std::stack< int > & stack_cols,
                       std::vector<bool> & visited_cols)
{
    int num_cols = matrix.size2();
    for ( int col = 0; col < num_cols; ++col )
    {
        if ( (matrix[row][col] != 0) && (visited_cols[col] == false) )
        {
            stack_cols.push(col);
            visited_cols[col] = true;
        }
    }
}

void add_incident_rows(int                 col,
                       MatrixT &           matrix,
                       std::stack< int > & stack_rows,
                       std::vector<bool> & visited_rows)
{
    int num_rows = matrix.size2();
    for ( int row = 0; row < num_rows; ++row )
    {
        if ( (matrix[row][col] != 0) && (visited_rows[row] == false) )
        {
            stack_rows.push(row);
            visited_rows[row] = true;
        }
    }
}

template< class MatrixT >
static Block find_block(int initial_row,
                        MatrixT & matrix,
                        std::vector<bool> & visited_rows,
                        std::vector<bool> & visited_cols)
{
    Component comp_rows;
    Component comp_cols;

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
            comp_cols.pus_back(col);
            stack_cols.pop();
        }
    }
    return std::make_pair(comp_rows, comp_cols);
}

template< class MatrixT >
Block BlockFinder::operator() (MatrixType & matrix)
{
    BlockPartition block_part;
    
    std::vector<bool> visited_cols;
    
    int num_rows = matrix.size1();
    std::vector<bool> visited_rows(num_rows, false);
    for ( int r = 0; r < num_rows; ++r )
    {
        if ( visited_rows[r] == false )
        {
            block_part.push_back(find_block(r));
        }
    }
    
    return block;
}

static int main ()
{
    BlockFinder block_finder;
    block_finder();
    return 0;
}
