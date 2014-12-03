#include "kappa.hpp"

void print_cohomology_class( const uint32_t g, const uint32_t m, const int32_t p, const VectorQ v )
{
    std::string filename_pref =
            "./cache/differentials_parallel/q_" +
            std::to_string(g) + "_" +
            std::to_string(m) + "_";
    MatrixQ image  = load_from_file_bz2< MatrixQ >( filename_pref + std::to_string(p-1) + "_base_changes", false );
    MatrixQ kernel = load_from_file_bz2< MatrixQ >( filename_pref + std::to_string(p) + "_triangular", false );
    MatrixQ::DiagonalType diagonal = load_from_file_bz2< MatrixQ::DiagonalType >( filename_pref + std::to_string(p-1) + "_diagonal", false );
    
//    std::cout << "Image:" << std::endl << image << std::endl << std::endl;
//    std::cout << "Kernel:" << std::endl << kernel << std::endl << std::endl;
//    std::cout << "Diagonal:" << std::endl << diagonal << std::endl << std::endl;
    std::cout << "Cohomology class: " << v.homology_class( kernel, image, diagonal ) << std::endl;
}

void add_cell( const MonoBasis& basis, VectorQ& v, const Q& alpha, const Tuple& cell )
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

void test_a()
{
    MonoBasis basis_0 = load_basis( 0, 1, 2 );
    VectorQ v(basis_0.size());
    
    Tuple a(2,1);
    a[1] = Transposition( 2, 1 );
    
    std::cout << a << std::endl;
    
    add_cell(basis_0, v, 1, a);
    
    std::cout << v << std::endl;
    
    print_cohomology_class(0, 1, 2, v); 
}

void test_aa()
{
    MonoBasis basis_0 = load_basis( 0, 2, 4 );
    VectorQ v(basis_0.size());
    
    Tuple aa(4,2);
    aa[1] = Transposition( 2, 1 );
    aa[2] = Transposition( 4, 3 );
    
    std::cout << aa << std::endl;
    
    add_cell(basis_0, v, 1, aa);
    
    std::cout << v << std::endl;
    
    print_cohomology_class(0, 2, 4, v); 
}

void test_c()
{
    MonoBasis basis_0 = load_basis( 1, 0, 4 );
    VectorQ v(basis_0.size());
    
    Tuple c(4,2);
    c[1] = Transposition( 3, 1 );
    c[2] = Transposition( 4, 2 );
    
    std::cout << c << std::endl;
    
    add_cell(basis_0, v, 1, c);
    
    std::cout << v << std::endl;
    
    print_cohomology_class(1, 0, 4, v);    
}

void test_d()
{
    MonoBasis basis_1 = load_basis( 1, 0, 3 );
    VectorQ v(basis_1.size());
    
    Tuple d(3,2);
    d[1] = Transposition( 2, 1 );
    d[2] = Transposition( 3, 1 );
    
    std::cout << d << std::endl;
    
    add_cell(basis_1, v, 1, d);
    
    std::cout << v << std::endl;
    
    print_cohomology_class(1, 0, 3, v);
}

void test_cc()
{
    MonoBasis basis_2 = load_basis( 2, 0, 8 );
    VectorQ v(basis_2.size());
    
    Tuple cc(8,4);
    cc[1] = Transposition( 3, 1 );
    cc[2] = Transposition( 4, 2 );
    cc[3] = Transposition( 7, 5 );
    cc[4] = Transposition( 8, 6 );
    
    std::cout << cc << std::endl;
    
    add_cell(basis_2, v, 1, cc);
    
    std::cout << v << std::endl;
    
    print_cohomology_class(2, 0, 8, v);
}

int main( int argc, char** argv )
{
    //test_a();
    test_aa();
    test_c();
    test_d();
    test_cc();
    return 0;
}
