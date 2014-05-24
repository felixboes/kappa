#include "matrix_field.hpp"

template< class CoefficientT >
MatrixField<CoefficientT>::MatrixField() : data(), num_rows(0), num_cols(0)
{
}

template< class CoefficientT >
MatrixField<CoefficientT>::MatrixField( size_t number_rows, size_t number_cols ) : data( number_rows*number_cols, CoefficientT(0) ), num_rows(number_rows), num_cols(number_cols)
{
}

template< class CoefficientT >
void MatrixField<CoefficientT>::row_operation( size_t row_1, size_t row_2, size_t col )
{
    CoefficientT lambda( -at(row_1,col) / at(row_2, col) );
    for( size_t j = col; j < num_cols; ++j )
    {
        CoefficientT & a = this->operator()( row_2, j );
        a = lambda * a + at( row_1, j );
    }
}

template< class CoefficientT >
 CoefficientT & MatrixField<CoefficientT>::operator()( size_t i, size_t j )
 {
     return data[ i*num_cols + j ];
 }

template< class CoefficientT >
void MatrixField<CoefficientT>::resize (size_t size1, size_t size2, bool)
{
    num_rows = size1;
    num_cols = size2;
    data.assign( size1 * size2, CoefficientT(0) );
    diagonal.clear();
}
 
template< class CoefficientT >
const CoefficientT & MatrixField<CoefficientT>::at(size_t i, size_t j) const
{
    return data.at( i*num_cols + j );
}

template< class CoefficientT >
size_t MatrixField<CoefficientT>::size1() const
{
    return num_rows;
}

template< class CoefficientT >
size_t MatrixField<CoefficientT>::size2() const
{
    return num_cols;
}

template< class CoefficientT >
void MatrixField<CoefficientT>::clear()
{
    for( auto & it : data )
    {
        it = CoefficientT(0);
    }
}

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const MatrixField<CoefficientT> & matrix )
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< class CoefficientT >
MatrixFieldCSS<CoefficientT>::MatrixFieldCSS() : data(), sec_data(), num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
{
}

template< class CoefficientT >
MatrixFieldCSS<CoefficientT>::MatrixFieldCSS( MatrixFieldCSSInitialization ini, size_t num_rows1, size_t num_cols1, size_t num_rows2, size_t num_cols2 ) :
    data(), sec_data(),
    num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0),
    row_operation_funct( &ThisType::row_operation_main_and_secondary )
{
    if( ini == only_main )
    {
        data = MatrixStorageType( num_rows1*num_cols1, CoefficientT(0) );
        num_rows = num_rows1;
        num_cols = num_cols1;
    }
    else if( ini == only_secondary )
    {
        sec_data = MatrixStorageType( num_rows1*num_cols1, CoefficientT(0) );
        sec_num_rows = num_rows1;
        sec_num_cols = num_cols1;
    }
    else if( ini == both )
    {
        data = MatrixStorageType( num_rows1*num_cols1, CoefficientT(0) );
        num_rows = num_rows1;
        num_cols = num_cols1;
        sec_data = MatrixStorageType( num_rows2*num_cols2, CoefficientT(0) );
        sec_num_rows = num_rows2;
        sec_num_cols = num_cols2;
    }
}

template< class CoefficientT >
MatrixFieldCSS<CoefficientT>::MatrixFieldCSS( size_t number_rows, size_t number_cols ) :
    data( number_rows*number_cols, CoefficientT(0) ),
    sec_data(),
    num_rows(number_rows),
    num_cols(number_cols),
    sec_num_rows(0),
    sec_num_cols(0)
{
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::row_operation( size_t row_1, size_t row_2, size_t col )
{
   (this->*row_operation_funct)(row_1, row_2, col);
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::row_operation_main_and_secondary( size_t row_1, size_t row_2, size_t col )
{
    CoefficientT lambda( -at(row_1,col) / at(row_2, col) );
    
    // Process main matrix
    for( size_t j = col; j < num_cols; ++j )
    {
        CoefficientT & a = this->operator()( row_2, j );
        a = lambda * a + at( row_1, j );
    }
    
    // Process secondary matrix
    for( size_t j = 0; j < sec_num_cols; ++j )
    {
        CoefficientT & a = sec( row_2, j );
        a = lambda * a + sec_at( row_1, j );
    }
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::row_operation_secondary( size_t row_1, size_t row_2, size_t col )
{
    CoefficientT lambda( -sec_at(row_1,col) / sec_at(row_2, col) );
    
    // Process secondary matrix
    for( size_t j = col; j < sec_num_cols; ++j )
    {
        CoefficientT & a = sec( row_2, j );
        a = lambda * a + sec_at( row_1, j );
    }
}

template< class CoefficientT >
void MatrixFieldCSS< CoefficientT >::define_row_operation( RowOperationType rt )
{
    if( rt == main_and_secondary )
    {
        row_operation_funct = &ThisType::row_operation_main_and_secondary;
    }
    else
    {
        row_operation_funct = &ThisType::row_operation_secondary;
    }
}

template< class CoefficientT >
CoefficientT & MatrixFieldCSS<CoefficientT>::operator()( size_t i, size_t j )
{
    return data[ i*num_cols + j ];
}

template< class CoefficientT >
CoefficientT & MatrixFieldCSS<CoefficientT>::sec( size_t i, size_t j )
{
    return sec_data[ i*sec_num_cols + j ];
}
 
template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::resize (size_t size1, size_t size2, bool)
{
    num_rows = size1;
    num_cols = size2;
    data.assign( size1 * size2, CoefficientT(0) );
    diagonal.clear();
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::sec_resize (size_t size1, size_t size2, bool)
{
    sec_num_rows = size1;
    sec_num_cols = size2;
    sec_data.assign( size1 * size2, CoefficientT(0) );
}

template< class CoefficientT >
const CoefficientT & MatrixFieldCSS<CoefficientT>::at(size_t i, size_t j) const
{
    return data.at( i*num_cols + j );
}

template< class CoefficientT >
const CoefficientT & MatrixFieldCSS<CoefficientT>::sec_at(size_t i, size_t j) const
{
    return sec_data.at( i*sec_num_cols + j );
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::size1() const
{
    return num_rows;
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::size2() const
{
    return num_cols;
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::sec_size1() const
{
    return sec_num_rows;
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::sec_size2() const
{
    return sec_num_cols;
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::clear()
{
    for( auto & it : data )
    {
        it = CoefficientT(0);
    }
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::sec_clear()
{
    for( auto & it : sec_data )
    {
        it = CoefficientT(0);
    }
}

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const MatrixFieldCSS<CoefficientT> & matrix )
{
    stream << "main " << matrix.num_rows << "x" << matrix.num_cols << ", secondary " << matrix.sec_num_rows << "x" << matrix.sec_num_cols << std::endl;
    stream << "main:" << std::endl;
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
