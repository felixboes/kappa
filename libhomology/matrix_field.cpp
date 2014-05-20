#include "matrix_field.hpp"

MatrixBool::MatrixBool() : data(), num_rows(0), num_cols(0)
{
}

MatrixBool::MatrixBool( size_t number_rows, size_t number_cols ) : data( number_rows, boost::dynamic_bitset<>(number_cols, 0) ), num_rows(number_rows), num_cols(number_cols)
{
}

void MatrixBool::row_operation( size_t row_1, size_t row_2, size_t )
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

 void MatrixBool::resize (size_t size1, size_t size2, bool)
 {
     for (size_t i = 0; i < data.size(); ++i)
     {
         data[i].reset();
     }

     data.resize(size1);
     for (size_t i = 0; i < size1; ++i)
     {
         data[i].resize(size2);
     }
     num_rows = size1;
     num_cols = size2;
 }
 
 bool MatrixBool::at(size_t i, size_t j) const
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

std::ostream& operator<< ( std::ostream& stream, const MatrixBool & matrix)
{
    for( size_t i = 0; i < matrix.num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.num_cols; )
        {
            stream << matrix.at(i,j);
            if( ++j < matrix.num_cols )
            {
                stream << ",";
            }
            else
            {
                stream << std::endl;
            }
        }
    }
    return stream;
}
