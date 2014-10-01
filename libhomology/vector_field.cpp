#include "vector_field.hpp"

VectorBool::VectorBool() : data(), dim(0)
{
}

VectorBool::VectorBool( size_t dimension ) : data( dimension, 0 ), dim(dimension)
{
}

bool VectorBool::operator()( const size_t i )
{
    return data[i];
}

void VectorBool::add_entry( const size_t i )
{
    data[i].flip();
}

void VectorBool::resize ( const size_t dimension, const bool )
{
    data.clear();   // this function removes all entries from the bitset. The result is a bitset of size 0.
    data.resize(dimension);
    dim = dimension;
}
 
bool VectorBool::at( const size_t i ) const
{
    return data[i];
}

size_t VectorBool::size() const
{
    return dim;
}

void VectorBool::clear()
{
    data.clear();   // this function removes all entries from the bitset. The result is a bitset of size 0.
    data.resize(dim);
}

std::ostream& operator<< ( std::ostream& stream, const VectorBool & vector )
{
    for( size_t i = 0; i < vector.dim - 1; ++i )
    {
        stream << std::setw(3) << vector.at(i) << ",";
    }
    stream << std::setw(3) << vector.at(vector.dim-1) << std::endl;
    return stream;
}
