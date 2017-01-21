#include "vector_field.hpp"
#include "vector_field_impl.ipp"

#include "homology.hpp"

/* Force template instantiation for used types */

#define force_template_instanciation( Coeff )\
    template class VectorField<Coeff>;\
    template std::ostream& operator<< ( std::ostream& stream, const VectorField<Coeff> & vector);\
    template void apply_base_changes_kernel( const MatrixField<Coeff>& m, VectorField<Coeff>& v );\
    template void apply_base_changes_image ( const MatrixField<Coeff>& m, VectorField<Coeff>& v );\
    template std::vector< VectorField<Coeff> > compute_base_of_kernel( const MatrixField<Coeff>& m );\
    template VectorField<Coeff> matrix_vector_product( const MatrixField<Coeff>& m, const VectorField<Coeff>& v );\
    bool matrix_vector_product_vanishes( const MatrixField<Coeff>& m, const VectorField<Coeff>& v );

force_template_instanciation(Q)
force_template_instanciation(Zm)

#undef force_template_instanciation

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

VectorBool& VectorBool::operator+=( const VectorBool& argument )
{
    data ^= argument.data;
    return *this;
}

VectorBool& VectorBool::operator-=( const VectorBool& argument )
{
    data ^= argument.data;
    return *this;
}

VectorBool& VectorBool::operator*=( const bool& argument )
{
    if( argument == false )
    {
        clear();
    }
    return *this;
}

void VectorBool::resize ( const size_t dimension, const bool )
{
    data.reset();
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
    data.reset();
}

std::ostream& operator<< ( std::ostream& stream, const VectorBool & vector )
{
    for( size_t i = 0; i < vector.dim - 1; ++i )
    {
        stream << std::setw(3) << (int32_t)vector.at(i) << ",";
    }
    stream << std::setw(3) << (int32_t)vector.at(vector.dim-1);
    return stream;
}

template<>
void apply_base_changes_kernel( const MatrixBool& m, VectorBool& v )
{
    size_t dim = v.size();
    const auto& diagonal = m.diagonal;
    
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
        for( size_t i = 0; i < dim; ++i )
        {
            // An entry does not encode a row operation if it is right of a diagonal element.
            if( diagonal_entry_occures_in_row[i] == true && diagonal_entry_col_row[i] <= diag_entry.second )
            {
            }
            else
            {
                if( m.at( i, diag_entry.second ) == true )
                {
                    v.add_entry(i);
                }
            }
        }
    }
}
