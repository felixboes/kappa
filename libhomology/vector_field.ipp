#include "vector_field.hpp"

template< class CoefficientT >
VectorField<CoefficientT>::VectorField() : data(), dim(0)
{
}

template< class CoefficientT >
VectorField<CoefficientT>::VectorField( const size_t dimension ) : data( dimension, CoefficientT(0) ), dim(dimension)
{
}

template< class CoefficientT >
CoefficientT & VectorField<CoefficientT>::operator()( const size_t i )
{
    return data[i];
}

template< class CoefficientT >
const CoefficientT & VectorField<CoefficientT>::at( const size_t i ) const
{
    return data.at(i);
}

template< class CoefficientT >
VectorField< CoefficientT >& VectorField< CoefficientT >::operator+=( const VectorField< CoefficientT >& argument )
{ 
    for( size_t i = 0; i < dim; ++i )
    {
        data[i] += argument.data[i];
    }
    return *this;
}

template< class CoefficientT >
VectorField< CoefficientT >& VectorField< CoefficientT >::operator-=( const VectorField< CoefficientT >& argument )
{
    for( size_t i = 0; i < dim; ++i )
    {
        data[i] -= argument.data[i];
    }
    return *this;
}

template< class CoefficientT >
VectorField< CoefficientT >& VectorField< CoefficientT >::operator*=( const CoefficientT& argument )
{
    for( size_t i = 0; i < dim; ++i )
    {
        data[i] *= argument;
    }
    return *this;
}

template< class CoefficientT >
void VectorField<CoefficientT>::resize( const size_t dimension, const bool )
{
    data.assign( dimension, CoefficientT(0) );
    dim = dimension;
}

template< class CoefficientT >
size_t VectorField<CoefficientT>::size() const
{
    return dim;
}

template< class CoefficientT >
void VectorField<CoefficientT>::clear()
{
    for( auto & it : data )
    {
        it = CoefficientT(0);
    }
}

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const VectorField<CoefficientT> & vector )
{
    stream << "[";
    for( size_t i = 0; i < vector.dim - 1; ++i )
    {
        stream << std::setw(3) << vector.at(i) << ",";
    }
    stream << std::setw(3) << vector.at(vector.dim - 1);
    return stream << "]";
}

template< class MatrixT, class VectorT >
void apply_base_changes( const MatrixT& m, VectorT& v )
{
    size_t dim = v.size();
    const auto& diagonal = m.diagonal;
    
    if( m.size1() == 0 || m.size2() == 0 )
    {
        return;
    }
    
    if( dim != m.size1() )
    {
        std::cout << "Error: The number of rows of the matrix is not equals the dimension of the vector." << std::endl;
        return;
    }
    if( diagonal.size() == 0 )
    {
        std::cout << "Error: The matrix seems to be not diagonalized. The diagonal of the matrix is empty." << std::endl;
        return;
    }

    // prepare fast access to the rows storing a diagonal entry.
    std::vector< bool >   diagonal_entry_occures_in_row (dim, false);
    std::vector< size_t > diagonal_entry_col_row (dim, 0);
    
    for( const auto& diag_entry : diagonal )
    {
        diagonal_entry_occures_in_row[diag_entry.first] = true;
        diagonal_entry_col_row[diag_entry.first] = diag_entry.second;
    }
    
    for( const auto& diag_entry : diagonal )
    {
        // compute new entries.
        const auto vector_entry = v.at( diag_entry.first );
        for( size_t i = 0; i < dim; ++i )
        {
            // An entry does not encode a row operation if it is right of a diagonal element.
            if( diagonal_entry_occures_in_row[i] == true && diagonal_entry_col_row[i] <= diag_entry.second )
            {
            }
            else
            {
                v(i) += m.at( i, diag_entry.second ) * vector_entry;
            }
        }
    }
}

template< class MatrixT, class VectorT >
VectorT matrix_vector_product( const MatrixT& m, const VectorT& v )
{
    size_t dim = v.size();
    VectorT res(dim);
    
    if( m.size1() == 0 || m.size2() == 0 )
    {
        return res;
    }
    if( dim != m.size2() )
    {
        std::cout << "Error: The number of rows of the matrix is not equals the dimension of the vector." << std::endl;
        return VectorT();
    }

    for( size_t i = 0; i < dim; ++i )
    {
        for( size_t j = 0; j < m.size2(); ++j )
        {
            res(i) += m.at(i,j) * v.at(j);
        }
    }
    return res;
}

template< class MatrixT, class VectorT >
bool matrix_vector_product_vanishes( const MatrixT& m, const VectorT& v )
{
    size_t dim = v.size();
    
    if( m.size1() == 0 || m.size2() == 0 )
    {
        return true;
    }
    if ( dim != m.size2() )
    {
        std::cout << "Error: The number of rows of the matrix is not equals the dimension of the vector." << std::endl;
        return false;
    }
    for( size_t i = 0; i < dim; ++i )
    {
        typename VectorT::CoefficientType res(0);
        for( size_t j = 0; j < m.size2(); ++j )
        {
            res += m.at(i,j) * v.at(j);
        }
        if( res != typename VectorT::CoefficientType(0) )
        {
            return false;
        }
    }
    return true;
}

