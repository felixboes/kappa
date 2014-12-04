#include "kappa.hpp"

template< class CoefficientT >
std::string filename_pref( uint32_t g, uint32_t m );

template<>
std::string filename_pref<Q>( uint32_t g, uint32_t m )
{
    return "./cache/differentials_parallel/q_" + std::to_string(g) + "_" + std::to_string(m) + "_";
}

template<>
std::string filename_pref<Zm>( uint32_t g, uint32_t m )
{
    return "./cache/differentials_parallel/s" + std::to_string( Zm::get_modulus() ) + "_" + std::to_string(g) + "_" + std::to_string(m) + "_";
}

template< class CoefficientT >
VectorField< CoefficientT > cohomology_class( const uint32_t g, const uint32_t m, const int32_t p, const VectorField< CoefficientT > v )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    std::cout << "Image:" << std::endl << image << std::endl << std::endl;
//    std::cout << "Kernel:" << std::endl << kernel << std::endl << std::endl;
//    std::cout << "Diagonal:" << std::endl << diagonal << std::endl << std::endl;
    std::cout << "The chain is " << ( matrix_vector_product_vanishes(load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false ), v) == true ? "indeed " : "NOT " ) << "a cycle." << std::endl;
    return v.homology_class( kernel, image );
}

template< class CoefficientT >
void add_cell( const MonoBasis& basis, VectorField< CoefficientT >& v, const CoefficientT& alpha, const Tuple& cell )
{
    int64_t id = basis.id_of(cell);
    
    if( id < 0 )
    {
        std::cout << "Error: The cell " << cell << " is not part of the basis." << std::endl;
        return;
    }
    
    v( (size_t)id ) += alpha;
}

MonoBasis load_basis( const uint32_t g, const uint32_t m, const int32_t p )
{
    std::string filename =
            "./cache/bases_parallel/" +
            std::to_string(g) + "_" +
            std::to_string(m) + "_" +
            std::to_string(p);
    
    return load_from_file_bz2< MonoBasis >( filename, false );
}

template< class CoefficientT >
void test_( const std::string& name, const uint32_t g, const uint32_t m, const uint32_t homological_p, const Tuple& cell )
{
    MonoBasis basis = load_basis( g, m, 4*g+2*m-homological_p );
    VectorField< CoefficientT > v(basis.size());
    add_cell<CoefficientT>(basis, v, 1, cell);
    
    std::cout << "Name                  = " << name << std::endl;
    std::cout << "Cell                  = " << cell << std::endl;
//    std::cout << "Vector in given basis = " << v << std::endl;
    std::cout << "Cohomology class      = " << cohomology_class( g, m, 4*g+2*m-homological_p, v ) << std::endl;
    std::cout << std::endl;
}

template< class CoefficientT >
void test_( const std::string& name, const uint32_t g, const uint32_t m, const uint32_t homological_p, const std::list<Tuple>& list)
{
    MonoBasis basis = load_basis( g, m, 4*g+2*m-homological_p );
    VectorField< CoefficientT > v(basis.size());
    
    std::cout << "Name                  = " << name << std::endl;
    for( const auto& it : list )
    {
        add_cell<CoefficientT>(basis, v, 1, it);
        std::cout << "Cell                  = " << it << std::endl;
    }
//    std::cout << "Vector in given basis = " << v << std::endl;
    std::cout << "Cohomology class      = " << cohomology_class( g, m, 4*g+2*m-homological_p, v ) << std::endl;
    std::cout << std::endl;
}

template< class CoefficientT >
void test_a()
{
    Tuple a(2,1);
    a[1] = Transposition( 2, 1 );
    test_<CoefficientT>("a", 0, 1, 0, a);
}

template< class CoefficientT >
void test_aa()
{
    Tuple aa(4,2);
    aa[1] = Transposition( 2, 1 );
    aa[2] = Transposition( 4, 3 );
    test_<CoefficientT>("a^2", 0, 2, 0, aa);
}

template< class CoefficientT >
void test_b()
{
    std::list< Tuple > list;
    
    Tuple cell(3,2);
    cell[1] = Transposition( 2, 1 );
    cell[2] = Transposition( 3, 2 );
    list.push_back(cell);
    
    cell[1] = Transposition( 3, 2 );
    cell[2] = Transposition( 3, 1 );
    list.push_back(cell);
    
    test_<CoefficientT>( "b", 0, 2, 1, list );
}

template< class CoefficientT >
void test_c()
{
    Tuple c(4,2);
    c[1] = Transposition( 3, 1 );
    c[2] = Transposition( 4, 2 );
    test_<CoefficientT>( "c", 1, 0, 0, c );
}

template< class CoefficientT >
void test_d()
{
    Tuple d(3,2);
    d[1] = Transposition( 2, 1 );
    d[2] = Transposition( 3, 1 );
    test_<CoefficientT>( "d", 1, 0, 1, d );
}

template< class CoefficientT >
void test_dd()
{
    Tuple dd(6,4);
    dd[1] = Transposition( 2, 1 );
    dd[2] = Transposition( 3, 1 );
    dd[3] = Transposition( 5, 4 );
    dd[4] = Transposition( 6, 4 );
    test_<CoefficientT>( "d^2", 2, 0, 2, dd );
}

template< class CoefficientT >
void test_aad()
{
    Tuple aad(7,4);
    aad[1] = Transposition( 2, 1 );
    aad[2] = Transposition( 3, 1 );
    aad[3] = Transposition( 5, 4 );
    aad[4] = Transposition( 7, 6 );
    test_<CoefficientT>( "aad", 1, 2, 1, aad );
}

template< class CoefficientT >
void test_bc()
{
    std::list< Tuple > list;
    
    Tuple cell(7,4);
    cell[1] = Transposition( 2, 1 );
    cell[2] = Transposition( 3, 2 );
    cell[3] = Transposition( 6, 4 );
    cell[4] = Transposition( 7, 5 );
    list.push_back(cell);
    
    cell[1] = Transposition( 3, 2 );
    cell[2] = Transposition( 3, 1 );
    cell[3] = Transposition( 6, 4 );
    cell[4] = Transposition( 7, 5 );
    list.push_back(cell);
    
    test_<CoefficientT>( "bc", 1, 2, 1, list );
}

template< class CoefficientT >
void test_cc()
{
    Tuple cc(8,4);
    cc[1] = Transposition( 3, 1 );
    cc[2] = Transposition( 4, 2 );
    cc[3] = Transposition( 7, 5 );
    cc[4] = Transposition( 8, 6 );
    test_<CoefficientT>( "c^2", 2, 0, 0, cc );
}

template< class CoefficientT >
void test_cd()
{
    Tuple cd(7,4);
    cd[1] = Transposition( 3, 1 );
    cd[2] = Transposition( 4, 2 );
    cd[3] = Transposition( 6, 5 );
    cd[4] = Transposition( 7, 5 );
    test_<CoefficientT>( "cd", 2, 0, 1, cd );
}


int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << "Program version: " << program_version_by_git << std::endl
              << "Date: " << current_date() << std::endl
              << std::endl;
    
    std::cout << "Rational computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    test_aa<Q>();
    test_b<Q>();
    test_c<Q>();
    test_d<Q>();
    test_dd<Q>();
    test_aad<Q>();
    test_bc<Q>();
    test_cc<Q>();
    test_cd<Q>();
    
    std::cout << "Mod 2 computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    Zm::set_modulus(2);
    test_aa<Zm>();
    test_b<Zm>();
    test_c<Zm>();
    test_d<Zm>();
    test_dd<Zm>();
    test_aad<Zm>();
    test_bc<Zm>();
    test_cc<Zm>();
    test_cd<Zm>();
    
    std::cout << "Mod 5 computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    Zm::set_modulus(5);
    test_bc<Zm>();
    test_dd<Zm>();
    test_cd<Zm>();
    
    return 0;
}
