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
    for( size_t i = 0; i < vector.dim - 1; ++i )
    {
        stream << std::setw(3) << vector.at(i) << ",";
    }
    stream << std::setw(3) << vector.at(vector.dim - 1) << std::endl;
    return stream;
}
