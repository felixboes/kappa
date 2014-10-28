#include "operationtester.hpp"

template< class MatrixComplex, class VectorT >
OperationTester< MatrixComplex, VectorT >::OperationTester( std::string coeff_prefix ) :
    coefficient_prefix( coeff_prefix )
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
            basis.insert( std::make_pair( idx, load_from_file_bz2< MonoBasis >( filename, print_status_messages ) ) );
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
            std::cout << "File '" << filename << ".bz2' does not exist. Cache generated?" << std::endl;
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
    forget_basis( MonoIndex( radial, genus, num_punctures, p ) ) ;
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

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT >::load_base_changes( const MonoIndex& idx, bool print_status_messages )
{
    std::string filename = 
        std::string("./cache/differentials_") + std::string( (std::get<0>(idx) == true ? "parallel" : "radial") ) + std::string("/") + 
        coefficient_prefix + "_" +
        std::to_string( std::get<1>(idx) ) + "_" +
        std::to_string( std::get<2>(idx) ) + "_" + 
        std::to_string( std::get<3>(idx) ) + "_base_changes";
    
    if( file_exists(filename + ".bz2") == true )
    {
        if( base_changes.count( idx ) == 0 )
        {
            base_changes.insert( std::make_pair( idx, load_from_file_bz2< MatrixType >( filename, print_status_messages ) ) );
        }
        else if( print_status_messages == true )
        {
            std::cout << "Base changes already exists." << std::endl;
        }
        return true;
    }
    else
    {
        if( print_status_messages == true )
        {
            std::cout << "File '" << filename << ".bz2' does not exist. Cache generated?" << std::endl;
        }
        return false;
    }
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT > :: load_base_changes( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages )
{
    return load_base_changes( MonoIndex(radial, genus, num_punctures, p), print_status_messages );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_base_changes( const MonoIndex& idx )
{
    base_changes.erase( idx );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_base_changes( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p )
{
    forget_base_changes( MonoIndex( radial, genus, num_punctures, p ) ) ;
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT >::load_triangular( const MonoIndex& idx, bool print_status_messages )
{
    std::string filename = 
        std::string("./cache/differentials_") + std::string( (std::get<0>(idx) == true ? "parallel" : "radial") ) + std::string("/") + 
        coefficient_prefix + "_" +
        std::to_string( std::get<1>(idx) ) + "_" +
        std::to_string( std::get<2>(idx) ) + "_" + 
        std::to_string( std::get<3>(idx) ) + "_triangular";
    
    if( file_exists(filename + ".bz2") == true )
    {
        if( triangular.count( idx ) == 0 )
        {
            triangular.insert( std::make_pair( idx, load_from_file_bz2< MatrixType >( filename, print_status_messages ) ) );
        }
        else if( print_status_messages == true )
        {
            std::cout << "Triangular already exists." << std::endl;
        }
        return true;
    }
    else
    {
        if( print_status_messages == true )
        {
            std::cout << "File '" << filename << ".bz2' does not exist. Cache generated?" << std::endl;
        }
        return false;
    }
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT > :: load_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages )
{
    return load_triangular( MonoIndex(radial, genus, num_punctures, p), print_status_messages );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_triangular( const MonoIndex& idx )
{
    triangular.erase( idx );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_triangular( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p )
{
    forget_triangular( MonoIndex( radial, genus, num_punctures, p ) ) ;
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT >::load_diagonal( const MonoIndex& idx, bool print_status_messages )
{
    std::string filename = 
        std::string("./cache/differentials_") + std::string( (std::get<0>(idx) == true ? "parallel" : "radial") ) + std::string("/") + 
        coefficient_prefix + "_" +
        std::to_string( std::get<1>(idx) ) + "_" +
        std::to_string( std::get<2>(idx) ) + "_" + 
        std::to_string( std::get<3>(idx) ) + "_diagonal";
    
    if( file_exists(filename + ".bz2") == true )
    {
        if( diagonal.count( idx ) == 0 )
        {
            diagonal.insert( std::make_pair( idx, load_from_file_bz2< DiagonalType >( filename, print_status_messages ) ) );
        }
        else if( print_status_messages == true )
        {
            std::cout << "Diagonal already exists." << std::endl;
        }
        return true;
    }
    else
    {
        if( print_status_messages == true )
        {
            std::cout << "File '" << filename << ".bz2' does not exist. Cache generated?" << std::endl;
        }
        return false;
    }
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT > :: load_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p, bool print_status_messages )
{
    return load_diagonal( MonoIndex(radial, genus, num_punctures, p), print_status_messages );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_diagonal( const MonoIndex& idx )
{
    diagonal.erase( idx );
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: forget_diagonal( bool radial, uint32_t genus, uint32_t num_punctures, int32_t p )
{
    forget_diagonal( MonoIndex( radial, genus, num_punctures, p ) ) ;
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: print_cache_status() const
{
    std::cout << std::endl
              << "---------------------- Cache Status ----------------------" << std::endl
              << "Memory Usage: " << current_memory_usage_in_mb() << " MB" << std::endl
              << std::endl
              << "Bases (" << basis.size() << "):" << std::endl;
    for( const auto& it : basis )
    {
        std::cout << "    " << it.first << std::endl;
    }
    
    std::cout << std::endl
              << "Base changes (" << base_changes.size() << "):" << std::endl;
    for( const auto& it : base_changes )
    {
        std::cout << "    " << it.first << std::endl;
    }
    
    std::cout << std::endl
              << "Triangular (" << triangular.size() << "):" << std::endl;
    for( const auto& it : triangular )
    {
        std::cout << "    " << it.first << std::endl;
    }
    
    std::cout << std::endl
              << "Diagonal (" << diagonal.size() << "):" << std::endl;
    for( const auto& it : diagonal )
    {
        std::cout << "    " << it.first << std::endl;
    }
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT > :: vector_is_valid( const MonoIndex& idx, const VectorType& v ) const
{
    if( basis.count( idx ) != 0 && basis.at( idx ).size() == v.size() )
    {
        return true;
    }
    else
    {
        return false;
    }
}

template< class MatrixComplex, class VectorT >
bool OperationTester< MatrixComplex, VectorT > :: vector_is_cycle( const MonoIndex& idx, const VectorType& v )
{
    if( triangular.count( idx ) == 0 )
    {
        std::cout << "Triangular Matrix " << idx << " is not yet loaded. Trying to load.";
        if( load_triangular( idx, false ) == false )
        {
            std::cout << " Failure." << std::endl;
            return false;
        }
        std::cout << " Success." << std::endl;
    }
    return matrix_vector_product_vanishes( triangular[idx], v );
}


