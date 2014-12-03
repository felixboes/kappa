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
    typedef typename MatrixField< CoefficientT >::DiagonalType DiagonalType;
    DiagonalType diagonal = load_from_file_bz2< DiagonalType >( filename_pref<CoefficientT>(g,m) + std::to_string(p-1) + "_diagonal", false );
//    std::cout << "Image:" << std::endl << image << std::endl << std::endl;
//    std::cout << "Kernel:" << std::endl << kernel << std::endl << std::endl;
//    std::cout << "Diagonal:" << std::endl << diagonal << std::endl << std::endl;
    return v.homology_class( kernel, image, diagonal );
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
void test_( const uint32_t g, const uint32_t m, const uint32_t homological_p, const Tuple& cell )
{
    MonoBasis basis = load_basis( g, m, 4*g+2*m-homological_p );
    VectorField< CoefficientT > v(basis.size());
    add_cell<CoefficientT>(basis, v, 1, cell);
    
    std::cout << "Cell                  = " << cell << std::endl;
    std::cout << "Vector in given basis = " << v << std::endl;
    std::cout << "Cohomology class      = " << cohomology_class( g, m, 4*g+2*m-homological_p, v ) << std::endl;
    std::cout << std::endl;
}

template< class CoefficientT >
void test_a()
{
    Tuple a(2,1);
    a[1] = Transposition( 2, 1 );
    test_<CoefficientT>(0, 1, 0, a);
}

template< class CoefficientT >
void test_aa()
{
    Tuple aa(4,2);
    aa[1] = Transposition( 2, 1 );
    aa[2] = Transposition( 4, 3 );
    test_<CoefficientT>(0, 2, 0, aa);
}

template< class CoefficientT >
void test_c()
{
    Tuple c(4,2);
    c[1] = Transposition( 3, 1 );
    c[2] = Transposition( 4, 2 );
    test_<CoefficientT>(1, 0, 0, c);
}

template< class CoefficientT >
void test_d()
{
    Tuple d(3,2);
    d[1] = Transposition( 2, 1 );
    d[2] = Transposition( 3, 1 );
    test_<CoefficientT>(1, 0, 1, d);
}

template< class CoefficientT >
void test_cc()
{
    Tuple cc(8,4);
    cc[1] = Transposition( 3, 1 );
    cc[2] = Transposition( 4, 2 );
    cc[3] = Transposition( 7, 5 );
    cc[4] = Transposition( 8, 6 );
    test_<CoefficientT>(2, 0, 0, cc);
}

int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << "Rational computations." << std::endl;
    test_aa<Q>();
    test_c<Q>();
    test_d<Q>();
    test_cc<Q>();
    
    std::cout << "Mod 2 computations." << std::endl;
    Zm::set_modulus(2);
    test_aa<Zm>();
    test_c<Zm>();
    test_d<Zm>();
    test_cc<Zm>();
    
    return 0;
}
