#include "matrix_field.hpp"

MatrixBool::MatrixBool() : diagonal(), data(), num_rows(0), num_cols(0)
{
}

MatrixBool::MatrixBool( size_t number_rows, size_t number_cols ) : diagonal(), data( number_rows, boost::dynamic_bitset<>(number_cols, 0) ), num_rows(number_rows), num_cols(number_cols)
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

 void MatrixBool::resize (size_t size1, size_t size2, bool )
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
     
     diagonal.clear();
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
        it.clear();
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

/////////////////////////////////////////////////////////////////////////////////////////////////

MatrixBoolCSS::MatrixBoolCSS() : diagonal(), data(), sec_data(), num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
{
}

MatrixBoolCSS::MatrixBoolCSS( const MatrixBoolCSSInitialization ini, const size_t num_rows1, const size_t num_cols1, const size_t num_rows2, const size_t num_cols2 ) :
    diagonal(), 
    data(), sec_data(),
    num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
{
    if( ini == only_main )
    {
        data = MatrixStorageType( num_rows1, boost::dynamic_bitset<>(num_cols1, 0) );
        num_rows = num_rows1;
        num_cols = num_cols1;
    }
    else if( ini == only_secondary )
    {
        sec_data = MatrixStorageType( num_rows1, boost::dynamic_bitset<>(num_cols1, 0) );
        sec_num_rows = num_rows1;
        sec_num_cols = num_cols1;
    }
    else if( ini == both )
    {
        data = MatrixStorageType( num_rows1, boost::dynamic_bitset<>(num_cols1, 0) );
        num_rows = num_rows1;
        num_cols = num_cols1;
        sec_data = MatrixStorageType( num_rows2, boost::dynamic_bitset<>(num_cols2, 0) );
        sec_num_rows = num_rows2;
        sec_num_cols = num_cols2;
    }
}

MatrixBoolCSS::MatrixBoolCSS( const size_t number_rows, const size_t number_cols ) :
    data( number_rows, boost::dynamic_bitset<>(number_cols, 0) ),
    sec_data(),
    num_rows(number_rows),
    num_cols(number_cols),
    sec_num_rows(0),
    sec_num_cols(0)
{
}

void MatrixBoolCSS::define_operations( const OperationType rt )
{
    if( rt == main_and_secondary )
    {
        row_operation_funct = &ThisType::row_operation_main_and_secondary;
        op_funct = &ThisType::main_op;
        at_funct = &ThisType::main_at;
        add_entry_funct = &ThisType::main_add_entry;
        size1_funct = &ThisType::main_size1;
        size2_funct = &ThisType::main_size2;
    }
    else
    {
        row_operation_funct = &ThisType::row_operation_secondary;
        op_funct = &ThisType::sec_op;
        at_funct = &ThisType::sec_at;
        add_entry_funct = &ThisType::sec_add_entry;
        size1_funct = &ThisType::sec_size1;
        size2_funct = &ThisType::sec_size2;
    }
}

bool MatrixBoolCSS::operator()( const size_t i, const size_t j )
{
   return (this->*op_funct)( i, j );
}


bool MatrixBoolCSS::main_op( const size_t i, const size_t j )
{
    return data[i][j];
}


bool MatrixBoolCSS::sec_op( const size_t i, const size_t j )
{
    return sec_data[i][j];
}


bool MatrixBoolCSS::at( const size_t i, const size_t j ) const
{
   return (this->*at_funct)( i, j );
}


bool MatrixBoolCSS::main_at(const size_t i, const size_t j) const
{
    return data[i][j];
}


bool MatrixBoolCSS::sec_at(const size_t i, const size_t j) const
{
    return sec_data[i][j];
}

void MatrixBoolCSS::row_operation( const size_t row_1, const size_t row_2, const size_t col )
{
   (this->*row_operation_funct)(row_1, row_2, col);
}

void MatrixBoolCSS::row_operation_main_and_secondary( const size_t row_1, const size_t row_2, const size_t )
{
    data[row_2] ^= data[row_1];
    sec_data[row_2] ^= sec_data[row_1];
}

void MatrixBoolCSS::row_operation_secondary( const size_t row_1, const size_t row_2, const size_t )
{
    sec_data[row_2] ^= sec_data[row_1];
}

void MatrixBoolCSS::add_entry( const size_t i, const size_t j )
{
    (this->*add_entry_funct)(i, j);
}

void MatrixBoolCSS::main_add_entry(const size_t i, const size_t j)
{
    data[i][j].flip();
}

void MatrixBoolCSS::sec_add_entry(const size_t i, const size_t j)
{
    sec_data[i][j].flip();
}

void MatrixBoolCSS::main_set( const size_t i, const size_t j, const bool val )
{
    data[i].set(j,val);
}

void MatrixBoolCSS::sec_set( const size_t i, const size_t j, const bool val )
{
    sec_data[i].set(j,val);
}

MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::main_row( const size_t i )
{
    return data[i];
}

const MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::main_row_at( const size_t i ) const
{
    return data.at(i);
}

MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::sec_row( const size_t i )
{
    return sec_data[i];
}

const MatrixBoolCSS::MatrixRowType& MatrixBoolCSS::sec_row_at( const size_t i ) const
{
    return sec_data.at(i);
}

void MatrixBoolCSS::resize(size_t size1, size_t size2, bool)
{
    for(size_t i = 0; i < data.size(); ++i)
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

void MatrixBoolCSS::sec_resize (const size_t size1, const size_t size2, const bool)
{
    for(size_t i = 0; i < sec_data.size(); ++i)
    {
        sec_data[i].reset();
    }
    
    sec_data.resize(size1);
    for (size_t i = 0; i < size1; ++i)
    {
        sec_data[i].resize(size2);
    }
    sec_num_rows = size1;
    sec_num_cols = size2;
}

size_t MatrixBoolCSS::size1() const
{
    return (this->*size1_funct)();
}

size_t MatrixBoolCSS::size2() const
{
    return (this->*size2_funct)();
}

size_t MatrixBoolCSS::main_size1() const
{
    return num_rows;
}

size_t MatrixBoolCSS::main_size2() const
{
    return num_cols;
}

size_t MatrixBoolCSS::sec_size1() const
{
    return sec_num_rows;
}

size_t MatrixBoolCSS::sec_size2() const
{
    return sec_num_cols;
}

void MatrixBoolCSS::clear()
{
    for( auto & it : data )
    {
        it.clear();
    }
}

void MatrixBoolCSS::sec_clear()
{
    for( auto & it : sec_data )
    {
        boost::dynamic_bitset<> & row = it;
        for ( size_t j = 0; j < sec_num_cols; ++j)
        {
            row[j] = 0;
        }
    }
}

std::ostream& operator<< ( std::ostream& stream, const MatrixBoolCSS & matrix )
{
    stream << "main " << matrix.num_rows << "x" << matrix.num_cols << ", secondary " << matrix.sec_num_rows << "x" << matrix.sec_num_cols << std::endl;
    stream << "main:" << std::endl;
    for( size_t i = 0; i < matrix.num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.num_cols; )
        {
            stream << matrix.main_at(i,j);
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
    
    stream << "secondary:" << std::endl;
    for( size_t i = 0; i < matrix.sec_num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.sec_num_cols; )
        {
            stream << matrix.sec_at(i,j);
            if( ++j < matrix.sec_num_cols )
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
