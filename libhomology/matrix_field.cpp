#include "matrix_field.hpp"

MatrixBool::MatrixBool() : data(), num_cols(0), num_rows(0)
{
}

MatrixBool::MatrixBool( size_t number_rows, size_t number_cols ) : data( number_rows, boost::dynamic_bitset<>(number_cols, 0) ), num_rows(number_rows), num_cols(number_cols)
{
}

void MatrixBool::row_operation( size_t row_1, size_t row_2, size_t col )
{
    data[row_2] ^= data[row_1];
}

 bool MatrixBool::operator()( size_t i, size_t j )
 {
     return data[i][j];
 }

 void MatrixBool::add_entry( size_t i, size_t j)
 {
        data[i][j].flip();
 }

 const bool MatrixBool::at(size_t i, size_t j) const
{
    return data[i][j];
}

size_t MatrixBool::size1() const
{
    return num_rows;
}

size_t MatrixBool::size2() const
{
    return num_cols;
}

void MatrixBool::clear()
{
    for( auto & it : data )
    {
        boost::dynamic_bitset<> & row = it;
        for ( size_t j = 0; j < num_cols; ++j)
        {
            row[j] = 0;
        }
    }
}

