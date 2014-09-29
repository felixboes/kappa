#include "matrix_field.hpp"

template< class CoefficientT >
MatrixField<CoefficientT>::MatrixField() : diagonal(), data(), num_rows(0), num_cols(0)
{
}

template< class CoefficientT >
MatrixField<CoefficientT>::MatrixField( const size_t number_rows, const size_t number_cols ) : diagonal(), data( number_rows*number_cols, CoefficientT(0) ), num_rows(number_rows), num_cols(number_cols)
{
}

template< class CoefficientT >
void MatrixField<CoefficientT>::row_operation( const size_t row_1, const size_t row_2, const size_t col )
{
    CoefficientT lambda( -at(row_2, col)/at(row_1,col)  );
    
    // Save lambda
    this->operator()( row_2, col ) = lambda;
    
    // Apply row operation on the remaining part
    for( size_t j = col+1; j < num_cols; ++j )
    {
        CoefficientT & a = this->operator()( row_2, j );
        a += lambda * at( row_1, j );
    }
    
}

template< class CoefficientT >
 CoefficientT & MatrixField<CoefficientT>::operator()( const size_t i, const size_t j )
 {
     return data[ i*num_cols + j ];
 }

template< class CoefficientT >
void MatrixField<CoefficientT>::resize (const size_t size1, const size_t size2, const bool)
{
    num_rows = size1;
    num_cols = size2;
    data.assign( size1 * size2, CoefficientT(0) );
    diagonal.clear();
}
 
template< class CoefficientT >
const CoefficientT & MatrixField<CoefficientT>::at(const size_t i, const size_t j) const
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
void MatrixField< CoefficientT > :: print_base_changes_in_short_form() const
{
    if( diagonal.size() == 0 )
    {
        std::cout << "The matrix seems to be not diagonalized: The diagonal of the matrix is empty." << std::endl;
        return;      
    }
    
    const size_t num_rows = size1();
    for( size_t i = 0; i < num_rows; ++i )
    {
        auto diag_entry = diagonal.cbegin();
        // print elements befor the diagonal entry.
        while( diag_entry != diagonal.cend() && diag_entry->first != i )
        {
            std::cout << std::setw(3) << this->at( i, diag_entry->second );
            ++diag_entry;
        }
        // print diagonal entry.
        if( diag_entry != diagonal.cend() )
        {
            std::cout << std::setw(3) << this->at( diag_entry->first, diag_entry->second );
            ++diag_entry;
        }
        // print elements after diagonal entry.
        while( diag_entry != diagonal.cend() )
        {
            std::cout << std::setw(3) << 0;
            ++diag_entry;
        }
        std::cout << std::endl;
    }
}

template< class CoefficientT >
void MatrixField< CoefficientT > :: print_triangular_shape() const
{
    if( diagonal.size() == 0 )
    {
        std::cout << "The matrix seems to be not diagonalized: The diagonal of the matrix is empty." << std::endl;
        return;      
    }
    
    std::cout << "Todo: Be clever." << std::endl;
}

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const MatrixField<CoefficientT> & matrix )
{
    // print diagonal if any.
    if( matrix.diagonal.size() > 0 )
    {
        std::cout << "Diagonal: ";
        for( const auto & diag_entry : matrix.diagonal )
        {
            std::cout << "(" << std::setw(2) << diag_entry.first << "," << std::setw(2) << diag_entry.second << ") ";
        }
        std::cout << std::endl;
    }
    
    // print matrix.
    for( size_t i = 0; i < matrix.num_rows; ++i )
    {
        for( size_t j = 0; j < matrix.num_cols; )
        {
            stream << std::setw(3) << matrix.at(i,j);
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
MatrixFieldCSS<CoefficientT>::MatrixFieldCSS() : diagonal(), data(), sec_data(), num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
{
}

template< class CoefficientT >
MatrixFieldCSS<CoefficientT>::MatrixFieldCSS( const MatrixFieldCSSInitialization ini, const size_t num_rows1, const size_t num_cols1, const size_t num_rows2, const size_t num_cols2 ) :
    diagonal(), 
    data(), sec_data(),
    num_rows(0), num_cols(0), sec_num_rows(0), sec_num_cols(0)
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
MatrixFieldCSS<CoefficientT>::MatrixFieldCSS( const size_t number_rows, const size_t number_cols ) :
    data( number_rows*number_cols, CoefficientT(0) ),
    sec_data(),
    num_rows(number_rows),
    num_cols(number_cols),
    sec_num_rows(0),
    sec_num_cols(0)
{
}

template< class CoefficientT >
void MatrixFieldCSS< CoefficientT >::define_operations( const OperationType rt )
{
    if( rt == main_and_secondary )
    {
        row_operation_funct = &ThisType::row_operation_main_and_secondary;
        op_funct = &ThisType::main_op;
        at_funct = &ThisType::main_at;
        size1_funct = &ThisType::main_size1;
        size2_funct = &ThisType::main_size2;
    }
    else
    {
        row_operation_funct = &ThisType::row_operation_secondary;
        op_funct = &ThisType::sec_op;
        at_funct = &ThisType::sec_at;
        size1_funct = &ThisType::sec_size1;
        size2_funct = &ThisType::sec_size2;
    }
}

template< class CoefficientT >
CoefficientT& MatrixFieldCSS<CoefficientT>::operator()( const size_t i, const size_t j )
{
   return (this->*op_funct)( i, j );
}

template< class CoefficientT >
CoefficientT & MatrixFieldCSS<CoefficientT>::main_op( const size_t i, const size_t j )
{
    return data[ i*num_cols + j ];
}

template< class CoefficientT >
CoefficientT & MatrixFieldCSS<CoefficientT>::sec_op( const size_t i, const size_t j )
{
    return sec_data[ i*sec_num_cols + j ];
}

template< class CoefficientT >
const CoefficientT& MatrixFieldCSS<CoefficientT>::at( const size_t i, const size_t j ) const
{
   return (this->*at_funct)( i, j );
}

template< class CoefficientT >
const CoefficientT & MatrixFieldCSS<CoefficientT>::main_at(const size_t i, const size_t j) const
{
    return data.at( i*num_cols + j );
}

template< class CoefficientT >
const CoefficientT & MatrixFieldCSS<CoefficientT>::sec_at(const size_t i, const size_t j) const
{
    return sec_data.at( i*sec_num_cols + j );
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::row_operation( const size_t row_1, const size_t row_2, const size_t col )
{
   (this->*row_operation_funct)(row_1, row_2, col);
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::row_operation_main_and_secondary( const size_t row_1, const size_t row_2, const size_t col )
{
    const CoefficientT lambda( -main_at(row_2, col)/main_at(row_1,col) );
    
    // Process main matrix
    for( size_t j = col; j < num_cols; ++j )
    {
        CoefficientT & a = main_op( row_2, j );
        a += lambda * main_at( row_1, j );
    }
    
    // Process secondary matrix
    for( size_t j = 0; j < sec_num_cols; ++j )
    {
        CoefficientT & a = sec_op( row_2, j );
        a += lambda * sec_at( row_1, j );
    }
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::row_operation_secondary( const size_t row_1, const size_t row_2, const size_t col )
{
    CoefficientT lambda( -sec_at(row_2, col) / sec_at(row_1,col) );
    
    // Process secondary matrix
    for( size_t j = col; j < sec_num_cols; ++j )
    {
        CoefficientT & a = sec_op( row_2, j );
        a += lambda * sec_at( row_1, j );
    }
}
 
template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::resize (const size_t size1, const size_t size2, const bool)
{
    num_rows = size1;
    num_cols = size2;
    data.assign( size1 * size2, CoefficientT(0) );
}

template< class CoefficientT >
void MatrixFieldCSS<CoefficientT>::sec_resize (const size_t size1, const size_t size2, const bool)
{
    sec_num_rows = size1;
    sec_num_cols = size2;
    sec_data.assign( size1 * size2, CoefficientT(0) );
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::size1() const
{
    return (this->*size1_funct)();
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::size2() const
{
    return (this->*size2_funct)();
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::main_size1() const
{
    return num_rows;
}

template< class CoefficientT >
size_t MatrixFieldCSS<CoefficientT>::main_size2() const
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
