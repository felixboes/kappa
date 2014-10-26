#include "operationtester.hpp"

template< class MatrixComplex, class VectorT >
OperationTester< MatrixComplex, VectorT >::OperationTester()
{
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT >::load_basis( const MonoIndex& idx, bool print_status_messages )
{
    std::string filename = 
        std::string("./cache/bases_") + std::string( (std::get<0>(idx) == true ? "parallel" : "radial") ) + std::string("/") + 
        std::to_string( std::get<1>(idx) ) + "_" +
        std::to_string( std::get<2>(idx) ) + "_" + 
        std::to_string( std::get<3>(idx) );
    
    if( file_exists(filename + ".bz2") == true )
    {
        if( basis.count( idx ) == 0 )
        {
            basis.insert( std::make_pair( idx, load_from_file_bz2< MonoBasis >( filename ) ) );
        }
        else if( print_status_messages == true )
        {
            std::cout << "Base already exists." << std::endl;
        }
        return true;
    }
    else
    {
        if( print_status_messages == true )
        {
            std::cout << "File '" << filename << "' does not exist. Cache generated?" << std::endl;
        }
        return false;
    }
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT > :: load_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages )
{
    return load_basis( MonoIndex(radial, genus, num_punctures, p), print_status_messages );
}

template< class MatrixComplex, class VectorT >
size_t OperationTester< MatrixComplex, VectorT > :: dim( const MonoIndex& idx ) const
{
    if( basis.count( idx ) == 0 )
    {
        return 0;
    }
    else
    {
        return basis.at( idx ).size();
    }
}

template< class MatrixComplex, class VectorT >
size_t OperationTester< MatrixComplex, VectorT > :: dim( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p ) const
{
    return dim( MonoIndex( radial, genus, num_punctures, p ) ) ;
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_basis( const MonoIndex& idx )
{
    basis.erase( idx );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_basis( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p )
{
    return forget_basis( MonoIndex( radial, genus, num_punctures, p ) ) ;
}

template< class MatrixComplex, class VectorT >
typename OperationTester< MatrixComplex, VectorT >::MonoIndex OperationTester< MatrixComplex, VectorT > :: product (
    const MonoIndex& idx_v,
    const MonoIndex& idx_w
)
{
    if( std::get<0>(idx_v) != std::get<0>(idx_w) )
    {
        return MonoIndex(false, 0, 0, 0);
    }
    
    return MonoIndex( std::get<0>(idx_v), std::get<1>(idx_v) + std::get<1>(idx_w), std::get<2>(idx_v) + std::get<2>(idx_w), std::get<3>(idx_v) + std::get<3>(idx_w) );
}
