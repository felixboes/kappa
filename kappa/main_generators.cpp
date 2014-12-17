#include "kappa.hpp"

MonoBasis load_basis( const uint32_t g, const uint32_t m, const int32_t p );

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
void cohomology_generators( const uint32_t g, const uint32_t m, const int32_t p )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MatrixField< CoefficientT > vanishing_test = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MonoBasis basis = load_basis( g, m, p );
    
    auto base = compute_base_of_kernel< MatrixField<CoefficientT>, VectorField<CoefficientT> >( kernel );
    for( const auto& v : base )
    {
        auto c = v.homology_class( kernel, image );
        if( c.is_zero() == false )
        {
            std::cout << "Chain = ";
            
            for( const auto& it : basis.basis )
            {
                if( v.at(it.id) != CoefficientT(0) )
                {
                    std::cout << std::setw(4) << v.at(it.id) << " * " << it << std::endl;
                }
            }
            
            std::cout << "The number of affiliated cells is " << v.number_non_vanishing_entries() << std::endl
                      << "The chain is " << ( matrix_vector_product_vanishes(vanishing_test, v) == true ? "indeed " : "NOT " ) << "a cycle." << std::endl
                      << "Its class it " << c << std::endl
                      << std::endl;
        }
    }
}

template< class CoefficientT >
void cohomology_generators_tex( const uint32_t g, const uint32_t m, const int32_t p )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MatrixField< CoefficientT > vanishing_test = load_from_file_bz2< MatrixField< CoefficientT > >( filename_pref<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MonoBasis basis = load_basis( g, m, p );
    
    
    std::cout << tex_preamble();
    auto base = compute_base_of_kernel< MatrixField<CoefficientT>, VectorField<CoefficientT> >( kernel );
    for( const auto& v : base )
    {
        auto c = v.homology_class( kernel, image );
        if( c.is_zero() == false )
        {
            std::cout << "The following chain represents the class " << c << " $\\in H_3(\\mathfrak M_{2,1}^0; \\mathbb K)$\\\\" << std::endl;
            
            for( const auto& it : basis.basis )
            {
                if( v.at(it.id) != CoefficientT(0) )
                {
                    //std::cout << "\\makebox[5ex][r]{$\\scriptstyle " << ( v.at(it.id) >= 0 ? "+" : "") << v.at(it.id) << "\\ $}";
                    std::cout << "\\makebox[5ex][r]{$\\scriptstyle " << v.at(it.id) << "\\ $}";
                    std::cout << tex_cell(it);
                }
            }
            
            std::cout << "\\clearpage" << std::endl;
        }
    }
    std::cout << tex_end();
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
void test_d_tex()
{
    Tuple d(3,2);
    d[1] = Transposition( 2, 1 );
    d[2] = Transposition( 3, 1 );
    std::cout << tex_cell(d);    
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
void test_dd_tex()
{
    Tuple dd(6,4);
    dd[1] = Transposition( 2, 1 );
    dd[2] = Transposition( 3, 1 );
    dd[3] = Transposition( 5, 4 );
    dd[4] = Transposition( 6, 4 );
    std::cout << tex_cell(dd);
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

template< class CoefficientT >
void test_z()
{
    std::list<Tuple> list;

    Tuple cell(5,4);
    cell[1] = Transposition( 4, 2 );
    cell[2] = Transposition( 4, 3 );
    cell[3] = Transposition( 4, 1 );
    cell[4] = Transposition( 5, 3 );
    list.push_back(cell);
    
    cell[1] = Transposition( 3, 2 );
    cell[2] = Transposition( 5, 1 );
    cell[3] = Transposition( 5, 4 );
    cell[4] = Transposition( 5, 2 );
    list.push_back(cell);
    
    cell[1] = Transposition( 3, 1 );
    cell[2] = Transposition( 4, 2 );
    cell[3] = Transposition( 4, 3 );
    cell[4] = Transposition( 5, 3 );
    list.push_back(cell);
    
    cell[1] = Transposition( 2, 1 );
    cell[2] = Transposition( 5, 3 );
    cell[3] = Transposition( 5, 1 );
    cell[4] = Transposition( 5, 4 );
    list.push_back(cell);
    
    test_<CoefficientT>( "z", 2, 0, 3, list );
}

template< class CoefficientT >
void test_z_candidates()
{
    std::vector<int> v_1 = {-1, -1,  1,  0,  1,  0,  0,  0,  1,  0,  1,  0, -1,  0,  0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  1,  0,  2,  0,  1,  0,  0,  0,  0,  2,  1,  1, -1,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0, -1,  0,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    std::vector<int> v_2 = { 0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  1,  0, -1, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  1,  0,  0,  2,  0,  1, -1,  1,  0,  0,  0,  0,  1,  0,  0,  0, -2,  0,  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  1,  1,  0, -1,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    std::vector<int> v_3 = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0};
    
    VectorField< CoefficientT > c_1(v_1.size());
    VectorField< CoefficientT > c_2(v_2.size());
    VectorField< CoefficientT > c_3(v_3.size());
    
    MonoBasis base = load_basis(2,0,5);
    
    for( size_t i = 0; i < v_1.size(); ++i )
    {
        c_1(i) = v_1[i];
        c_2(i) = v_2[i];
        c_3(i) = v_3[i];
    }

    std::cout << "Print cells of the first vector. There are " << c_1.number_non_vanishing_entries() << " cells involved." << std::endl;
    for( const auto& it : base.basis )
    {
        if( c_1.at(it.id) != CoefficientT(0) )
        {
            std::cout << std::setw(4) << c_1.at(it.id) << " * " << it << std::endl;
        }
    }
    
    std::cout << "Print cells of the second vector. There are " << c_2.number_non_vanishing_entries() << " cells involved." << std::endl;
    for( const auto& it : base.basis )
    {
        if( c_2.at(it.id) != CoefficientT(0) )
        {
            std::cout << std::setw(4) << c_2.at(it.id) << " * " << it << std::endl;
        }
    }
    
    std::cout << "Print cells of the second vector. There are " << c_2.number_non_vanishing_entries() << " cells involved." << std::endl;
    for( const auto& it : base.basis )
    {
        if( c_2.at(it.id) != CoefficientT(0) )
        {
            std::cout << std::setw(4) << c_2.at(it.id) << " * " << it << std::endl;
        }
    }
    
    std::cout << "Print cells of the third vector. There are " << c_3.number_non_vanishing_entries() << " cells involved." << std::endl;
    for( const auto& it : base.basis )
    {
        if( c_3.at(it.id) != CoefficientT(0) )
        {
            std::cout << std::setw(4) << c_3.at(it.id) << " * " << it << std::endl;
        }
    }
}

template< class CoefficientT >
void test_z_candidates_tex()
{
    std::vector<int> v_1 = {-1, -1,  1,  0,  1,  0,  0,  0,  1,  0,  1,  0, -1,  0,  0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  1,  0,  2,  0,  1,  0,  0,  0,  0,  2,  1,  1, -1,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0, -1,  0,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    std::vector<int> v_2 = { 0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  1,  0, -1, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  1,  0,  0,  2,  0,  1, -1,  1,  0,  0,  0,  0,  1,  0,  0,  0, -2,  0,  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  1,  1,  0, -1,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
    
    VectorField< CoefficientT > c_1(v_1.size());
    VectorField< CoefficientT > c_2(v_2.size());
    
    MonoBasis base = load_basis(2,0,5);
    
    for( size_t i = 0; i < v_1.size(); ++i )
    {
        c_1(i) = v_1[i];
        c_2(i) = v_2[i];
    }
    
    std::cout << tex_preamble();
    std::cout << "Print cells of the first vector. There are " << c_1.number_non_vanishing_entries() << " cells involved." << std::endl << std::endl;
    for( const auto& it : base.basis )
    {
        if( c_1.at(it.id) != CoefficientT(0) )
        {
            std::cout << tex_cell(it);
        }
    }
    
    std::cout << std::endl;
    
    std::cout << "Print cells of the second vector. There are " << c_2.number_non_vanishing_entries() << " cells involved." << std::endl << std::endl;
    for( const auto& it : base.basis )
    {
        if( c_2.at(it.id) != CoefficientT(0) )
        {
            std::cout << tex_cell(it);
        }
    }
    
    std::cout << tex_end();
}

int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << "Program version: " << program_version_by_git << std::endl
              << "GMP Version:     " << gmp_version << std::endl
              << "Boost Version:   " << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl
              << "Date:            " << current_date() << std::endl
              << std::endl;
    
//    std::cout << "Rational computations." << std::endl;
//    std::cout << "--------------------------------" << std::endl;
//    test_aa<Q>();
//    test_b<Q>();
//    test_c<Q>();
//    test_d<Q>();
//    test_dd<Q>();
//    test_aad<Q>();
//    test_bc<Q>();
//    test_cc<Q>();
//    test_cd<Q>();
//    test_z<Q>();
    
//    std::cout << "Mod 2 computations." << std::endl;
//    std::cout << "--------------------------------" << std::endl;
//    Zm::set_modulus(2);
//    test_aa<Zm>();
//    test_b<Zm>();
//    test_c<Zm>();
//    test_d<Zm>();
//    test_dd<Zm>();
//    test_aad<Zm>();
//    test_bc<Zm>();
//    test_cc<Zm>();
//    test_cd<Zm>();
//    test_z<Zm>();
    
//    std::cout << "Mod 5 computations." << std::endl;
//    std::cout << "--------------------------------" << std::endl;
//    Zm::set_modulus(5);
//    test_bc<Zm>();
//    test_dd<Zm>();
//    test_cd<Zm>();

//    cohomology_generators<Q>( 0, 2, 3);
//    cohomology_generators<Q>( 1, 0, 4);
//    cohomology_generators<Q>( 1, 0, 3);
    //cohomology_generators<Q>( 2, 0, 6);
    //cohomology_generators_tex<Q>( 2, 0, 5);

    //test_z_candidates<Q>();
    //test_z_candidates_tex<Q>();
    
    std::cout << tex_preamble();
    test_d_tex<Q>();
    test_dd_tex<Q>();
    std::cout << tex_end();
    
    return 0;
}
