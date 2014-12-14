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

template< class MatrixComplex, class VectorT >
typename OperationTester< MatrixComplex, VectorT >::VectorType OperationTester< MatrixComplex, VectorT > :: vector_homology_class( const MonoIndex& idx, const VectorType& v )
{
    size_t dim = v.size();
    
    if( vector_is_valid(idx, v) == false )
    {
        std::cout << "Vector " << v << " is not valid." << std::endl;
        return VectorT();
    }
    if( base_changes.count( idx ) == 0 )
    {
        std::cout << "Bases change " << idx << " is not yet loaded. Trying to load.";
        if( load_base_changes( idx, false ) == false )
        {
            std::cout << " Failure." << std::endl;
            return VectorT();
        }
        std::cout << " Success." << std::endl;
    }
    
    if( diagonal.count( idx ) == 0 )
    {
        std::cout << "Diagonal " << idx << " is not yet loaded. Trying to load.";
        if( load_diagonal( idx, false ) == false )
        {
            std::cout << " Failure." << std::endl;
            return VectorT();
        }
        std::cout << " Success." << std::endl;
    }
    
    // apply base changes.
    VectorT intermediate_homology_class = v;
    apply_base_changes( base_changes[idx], intermediate_homology_class );
    
    // compute betti number.
    size_t betti_number = dim;
    std::vector< bool > diagonal_entry_occures_in_row ( dim, false );
    for( const auto& diag_entry : diagonal[idx] )
    {
        diagonal_entry_occures_in_row[diag_entry.first] = true;
        --betti_number;
    }
    
    // project onto homology modules.
    VectorT homology_class(betti_number);
    size_t current_dim = 0;
    for( size_t i = 0; i < dim; ++i )
    {
        if( diagonal_entry_occures_in_row.at(i) == false )
        {
            homology_class(current_dim) = intermediate_homology_class(i);
            ++current_dim;
        }
    }
    return homology_class;
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
typename OperationTester< MatrixComplex, VectorT >::VectorType OperationTester< MatrixComplex, VectorT > :: product(
    const MonoIndex& idx_v,
    const VectorType& v,
    const MonoIndex& idx_w,
    const VectorType& w )
{
    if( vector_is_valid(idx_v, v) == false || vector_is_valid(idx_w, w) == false )
    {
        std::cout << "At least one of the given vectors is not valid." << std::endl;
        return VectorType();
    }
    MonoIndex idx_prod = product( idx_v, idx_w );
    
    if( basis.count(idx_prod) == 0 )
    {
        std::cout << "Basis " << idx_prod << " is not yet loaded. Trying to load.";
        if( load_base_changes( idx_prod, false ) == false )
        {
            std::cout << " Failure." << std::endl;
            return VectorType();
        }
        std::cout << " Success." << std::endl;
    }
    
    size_t dim_prod = basis.at(idx_prod).size();
    VectorType vect_prod(dim_prod);
    
    for( const auto& tuple_v : basis[idx_v].basis )
    {
        for( const auto& tuple_w : basis[idx_w].basis )
        {
            size_t j = basis[idx_prod].id_of( tuple_v * tuple_w );
            vect_prod(j) = v.at(tuple_v.id) * w.at(tuple_w.id);
        }
    }
    
    return vect_prod;
}

template< class MatrixComplex, class VectorT >
typename OperationTester< MatrixComplex, VectorT >::VectorType OperationTester< MatrixComplex, VectorT > :: Q( const MonoIndex& idx, const VectorType& v )
{
    MonoIndex idx_res( std::get<0>(idx), 2 * std::get<1>(idx), 2 * std::get<2>(idx), 2 * std::get<3> - 1 );
    load_basis(idx_res, false);
    load_basis(idx, false);
    size_t dim_res = dim(idx_res);
    VectorType res( dim_res );
    const CoefficientType zero(0);
    const auto basis_v( basis.at( idx ).basis );
    
    for( const auto& basis_element : basis_v )
    {
        if( v( basis_v( basis_element ) ) != zero )
        {
            compute_and_add_Q( basis_v.at(c), basis_element, v);
        }
    }
    
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: compute_and_add_Q( const CoefficientType& c, const Tuple& t, VectorType v)
{
    
}

template< class MatrixComplex, class VectorT >
void OperationTester< MatrixComplex, VectorT > :: compute_and_add_kappa_dual( const CoefficientType& c, const Tuple& t, VectorType v)
{
    
}
