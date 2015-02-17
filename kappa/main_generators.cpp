#include "kappa.hpp"

template< class CoefficientT >
VectorField< CoefficientT > cohomology_class( const uint32_t g, const uint32_t m, const int32_t p, const VectorField< CoefficientT > v )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    std::cout << "Image:" << std::endl << image << std::endl << std::endl;
//    std::cout << "Kernel:" << std::endl << kernel << std::endl << std::endl;
//    std::cout << "Diagonal:" << std::endl << diagonal << std::endl << std::endl;
    std::cout << "The chain is " << ( matrix_vector_product_vanishes(load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false ), v) == true ? "indeed " : "NOT " ) << "a cycle." << std::endl;
    return v.homology_class( kernel, image );
}

template< class CoefficientT >
void cohomology_generators( const uint32_t g, const uint32_t m, const int32_t p )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MatrixField< CoefficientT > vanishing_test = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MonoBasis basis = load_parallel_mono_basis( g, m, p );
    
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
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MatrixField< CoefficientT > vanishing_test = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
    MonoBasis basis = load_parallel_mono_basis( g, m, p );
    
    
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

template< class CoefficientT >
void test_( const std::string& name, const uint32_t g, const uint32_t m, const uint32_t homological_p, const std::list<Tuple>& list)
{
    MonoBasis basis = load_parallel_mono_basis( g, m, 4*g+2*m-homological_p );
    VectorField< CoefficientT > v(basis.size());
    
    std::cout << "Name                  = " << name << std::endl;
    for( const auto& cell : list )
    {
        v += kappa_dual< VectorField< CoefficientT > >( 1, cell, basis );
//        std::cout << "Cell                  = " << cell << std::endl;
//        std::cout << "Kappa Dual            = " << kappa_dual< VectorField< CoefficientT > >( 1, cell, basis ) << std::endl;
    }
//    std::cout << "Vector in given basis = " << v << std::endl;
    std::cout << "Cohomology class      = " << cohomology_class( g, m, 4*g+2*m-homological_p, v ) << std::endl;
    std::cout << std::endl;
}

template< class CoefficientT >
void test_( const std::string& name, const uint32_t g, const uint32_t m, const uint32_t homological_p, const Tuple& cell )
{
    std::list< Tuple > list;
    list.push_back(cell);
    
    test_< CoefficientT >(name, g, m, homological_p, list);
}

template< class CoefficientT >
void test_a()
{
    test_<CoefficientT>("a", 0, 1, 0, create_cell(1, 2, 1) );
}

template< class CoefficientT >
void test_aa()
{
    test_<CoefficientT>("a^2", 0, 2, 0, create_cell(2, 4, 3, 2, 1) );
    test_<CoefficientT>("a^2", 0, 2, 0, create_cell(2, 4, 1, 3, 2) );
}

template< class CoefficientT >
void test_b()
{
    std::list< Tuple > list;
    list.push_back( create_cell(2, 3, 2, 2, 1) );
    list.push_back( create_cell(2, 3, 1, 3, 2) );
    
    test_<CoefficientT>( "b", 0, 2, 1, list );
}

template< class CoefficientT >
void test_ab()
{
    std::list< Tuple > list;
    
    list.push_back( create_cell(3, 5, 4, 3, 2, 2, 1) );
    list.push_back( create_cell(3, 5, 4, 3, 1, 3, 2) );
    
    test_<CoefficientT>( "ab", 0, 3, 1, list );
}

template< class CoefficientT >
void test_bb()
{
    std::list< Tuple > list;
    std::list< Tuple > helper_list;
    
    helper_list.push_back( create_cell(2, 3, 2, 2, 1) );
    helper_list.push_back( create_cell(2, 3, 1, 3, 2) );
    
    list.push_back( helper_list.front() * helper_list.front() );
    list.push_back( helper_list.front() * helper_list.back() );
    list.push_back( helper_list.back()  * helper_list.front() );
    list.push_back( helper_list.back()  * helper_list.back() );
   
    test_<CoefficientT>( "bb", 0, 4, 2, list );
}

template< class CoefficientT >
void test_c()
{
    test_<CoefficientT>( "c", 1, 0, 0, create_cell(2, 4, 2, 3, 1) );
}

template< class CoefficientT >
void test_d()
{
    test_<CoefficientT>( "d", 1, 0, 1, create_cell(2, 3, 1, 2, 1) );
}

void test_d_tex()
{
    std::cout << tex_cell( create_cell(2, 3, 1, 2, 1) );    
}

template< class CoefficientT >
void test_dd()
{
    Tuple d = create_cell(2, 3, 1, 2, 1);
    test_<CoefficientT>( "dd", 2, 0, 2, d*d );
}

void test_dd_tex()
{
    Tuple d = create_cell(2, 3, 1, 2, 1);
    std::cout << tex_cell(d*d);
}

template< class CoefficientT >
void test_aad()
{
    Tuple aad = create_cell(4, 7, 6, 5, 4, 3, 1, 2, 1);
    test_<CoefficientT>( "aad", 1, 2, 1, aad );
}

template< class CoefficientT >
void test_bc()
{
    std::list< Tuple > list;
    
    list.push_back( create_cell(4, 7, 5, 6, 4, 3, 2, 2, 1) );
    list.push_back( create_cell(4, 7, 5, 6, 4, 3, 1, 3, 2) );
    
    test_<CoefficientT>( "bc", 1, 2, 1, list );
}

template< class CoefficientT >
void test_ccb()
{
    std::list< Tuple > list;
    
    Tuple c = create_cell(2, 4, 2, 3, 1);
    
    list.push_back( c * c * create_cell(2, 3, 2, 2, 1) );
    list.push_back( c * c * create_cell(2, 3, 1, 3, 2) );
    
    test_<CoefficientT>( "ccb", 2, 2, 1, list );
}

template< class CoefficientT >
void test_cc()
{
    Tuple c = create_cell(2, 4, 2, 3, 1);
    test_<CoefficientT>( "c^2", 2, 0, 0, c*c );
}

template< class CoefficientT >
void test_cd()
{
    test_<CoefficientT>( "cd", 2, 0, 1, create_cell(2, 4, 2, 3, 1) * create_cell(2, 3, 1, 2, 1) );
}

template< class CoefficientT >
void test_e()
{
    test_<CoefficientT>( "e", 1, 1, 2, create_cell(3, 3, 1, 4, 3, 2, 1)  );
}

template< class CoefficientT >
void test_ce()
{
    test_<CoefficientT>( "ce", 2, 1, 2, create_cell(2, 4, 2, 3, 1) * create_cell(3, 3, 1, 4, 3, 2, 1)  );
}

template< class CoefficientT >
void test_de()
{
    test_<CoefficientT>( "de", 2, 1, 3, create_cell(2, 3, 1, 2, 1) * create_cell(3, 3, 1, 4, 3, 2, 1) );
}

template< class CoefficientT >
void test_Qc()
{
    std::list< Tuple > list;
    
    list.push_back( create_cell(4, 7, 5, 6, 4, 4, 2, 3, 1) );
    list.push_back( create_cell(4, 7, 5, 6, 1, 4, 2, 3, 1) );
    list.push_back( create_cell(4, 7, 5, 6, 2, 4, 2, 3, 1) );
    list.push_back( create_cell(4, 7, 5, 6, 3, 4, 2, 3, 1) );
    
    list.push_back( create_cell(4, 7, 5, 6, 1, 5, 3, 4, 2) );
    list.push_back( create_cell(4, 7, 2, 6, 1, 5, 3, 4, 2) );
    list.push_back( create_cell(4, 7, 3, 6, 1, 5, 3, 4, 2) );
    list.push_back( create_cell(4, 7, 4, 6, 1, 5, 3, 4, 2) );
    
    list.push_back( create_cell(4, 7, 2, 6, 1, 6, 4, 5, 3) );
    list.push_back( create_cell(4, 7, 2, 3, 1, 6, 4, 5, 3) );
    list.push_back( create_cell(4, 7, 2, 4, 1, 6, 4, 5, 3) );
    list.push_back( create_cell(4, 7, 2, 5, 1, 6, 4, 5, 3) );
    
    list.push_back( create_cell(4, 7, 2, 3, 1, 7, 5, 6, 4) );
    list.push_back( create_cell(4, 4, 2, 3, 1, 7, 5, 6, 4) );
    list.push_back( create_cell(4, 5, 2, 3, 1, 7, 5, 6, 4) );
    list.push_back( create_cell(4, 6, 2, 3, 1, 7, 5, 6, 4) );
    
    test_<CoefficientT>( "Qc", 2, 0, 1, list );
}

std::list< Tuple > create_Qd()
{
    std::list<Tuple> list;
    list.push_back( create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 3, 4, 2, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 1, 4, 2, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 2, 4, 2, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 1, 4, 1, 4, 2, 3, 2 ) );
    list.push_back( create_cell( 4, 5, 2, 4, 1, 4, 2, 3, 2 ) );
    list.push_back( create_cell( 4, 5, 3, 4, 1, 4, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 5, 1, 2, 1, 4, 2, 3, 2 ) );
    list.push_back( create_cell( 4, 5, 3, 2, 1, 4, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 5, 1, 2, 1, 5, 3, 4, 3 ) );
    list.push_back( create_cell( 4, 3, 1, 2, 1, 5, 3, 4, 3 ) );
    list.push_back( create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
    
    return list;
}

template< class CoefficientT >
void test_Qd()
{
    test_<CoefficientT>( "Q(d)", 2, 0, 3, create_Qd() );
}

void test_Qd_tex()
{
    std::cout << tex_cell( create_Qd() );
}

template< class CoefficientT >
void test_Te()
{
    std::list< Tuple > list;
    list.push_back( create_cell(4, 5, 3, 3, 1, 4, 3, 2, 1) );
    list.push_back( create_cell(4, 5, 1, 3, 1, 4, 3, 2, 1) );
    
    test_<CoefficientT>( "T(e)", 2, 0, 3, list );
}

std::list< Tuple > create_z_1()
{
    std::list<Tuple> list;
    list.push_back( create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 3, 4, 2, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 1, 4, 2, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 2, 4, 2, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 1, 4, 1, 4, 2, 3, 2 ) );
    list.push_back( create_cell( 4, 5, 2, 4, 1, 4, 2, 3, 2 ) );
    list.push_back( create_cell( 4, 5, 3, 4, 1, 4, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 5, 1, 2, 1, 4, 2, 3, 2 ) );
    list.push_back( create_cell( 4, 5, 3, 2, 1, 4, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 4, 2, 3, 2, 5, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 4, 1, 3, 1, 5, 2, 3, 2 ) );
    
    list.push_back( create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
    
    return list;
}

template< class CoefficientT >
void test_z_1()
{
    test_<CoefficientT>( "z", 2, 0, 3, create_z_1() );
}

std::list< Tuple > create_z_2()
{
    std::list<Tuple> list;
    list.push_back( create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 1, 4, 3, 4, 1, 2, 1 ) );
    list.push_back( create_cell( 4, 5, 3, 4, 3, 4, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 5, 3, 3, 1, 4, 1, 2, 1 ) );
    
    list.push_back( create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
    
    return list;
}

template< class CoefficientT >
void test_z_2()
{
    test_<CoefficientT>( "z", 2, 0, 3, create_z_2() );
}

template< class CoefficientT >
void test_tilde_d()
{
    std::list< Tuple > list;
    list.push_back( create_cell(4, 5, 4, 5, 4, 3, 1, 2, 1) );
    list.push_back( create_cell(4, 5, 3, 4, 3, 2, 1, 2, 1) );
    list.push_back( create_cell(4, 4, 3, 5, 4, 2, 1, 2, 1) );
    list.push_back( create_cell(4, 5, 4, 5, 4, 2, 1, 3, 2) );
    test_<CoefficientT>( "tilde d", 2, 0, 3, list );
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
    
    MonoBasis base = load_parallel_mono_basis(2,0,5);
    
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
    
    MonoBasis base = load_parallel_mono_basis(2,0,5);
    
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

template< class CoefficientT >
void test_Q()
{
    Tuple a(2,1);
    a[1] = Transposition(2,1);
    
    const auto slits = a.slits();
    for( const auto& it : slits )
    {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    
    const auto shuffle_pos = a.shuffle_positions();
    for( const auto& it : shuffle_pos )
    {
        std::cout << it << " ";
    }
    std::cout << std::endl;
}


int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << kappa_version( argc, argv ) << std::endl;
    
    std::cout << create_cochain<Q>( Generator::a ) << std::endl;
    std::cout << create_cochain<Q>( Generator::b ) << std::endl;
    std::cout << create_cochain<Q>( Generator::c ) << std::endl;
    std::cout << create_cochain<Q>( Generator::d ) << std::endl;
    std::cout << create_cochain<Q>( Generator::e ) << std::endl;
    
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
//    test_Qd<Q>();
//    test_Te<Q>();
//    test_z_1<Q>();
//    test_z_2<Q>();
    
//    std::cout << "Mod 2 computations." << std::endl;
//    std::cout << "--------------------------------" << std::endl;
//    Zm::set_modulus(2);
//    test_aa<Zm>();
//    test_b<Zm>();
//    test_ab<Zm>();
//    test_bb<Zm>();
//    test_c<Zm>();
//    test_d<Zm>();
//    test_dd<Zm>();
//    test_aad<Zm>();
//    test_bc<Zm>();
//    test_cc<Zm>();
//    test_cd<Zm>();
//    test_Qd<Zm>();
//    test_Te<Zm>();
//    test_z_1<Zm>();
//    test_z_2<Zm>();
//    test_tilde_d<Zm>();
//    test_e<Zm>();
//    test_de<Zm>();
//    test_ccb<Zm>();
//    test_Qc<Zm>();
    
    //cohomology_generators<Zm>(2, 0, 2*2*2-3);
    //cohomology_generators<Zm>(2, 0, 4*2 + 2*0 - 2);
    
//    std::cout << "Mod 5 computations." << std::endl;
//    std::cout << "--------------------------------" << std::endl;
//    Zm::set_modulus(5);
//    test_bc<Zm>();
//    test_dd<Zm>();
//    test_cd<Zm>();
//    test_Qd<Zm>();
//    test_Te<Zm>();
//    test_Qc<Zm>();
    
//    cohomology_generators<Q>( 0, 2, 3);
//    cohomology_generators<Q>( 1, 0, 4);
//    cohomology_generators<Q>( 1, 0, 3);
//    cohomology_generators<Q>( 2, 1, 4);

//    test_z_candidates<Q>();
//    test_z_candidates_tex<Q>();
    
//    std::cout << tex_preamble();
//    test_d_tex<Q>();
//    test_dd_tex<Q>();
//    test_Qd_tex();
//    std::cout << tex_end();

    //test_Q<Q>();
   
//    const auto res1 = list_set_partitions(4);
    
//    for( const auto& it : res1 )
//    {
//        for( const auto& inner : it )
//        {
//            for( const auto& inner2 : inner )
//            {
//                std::cout << inner2 << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    
//    const auto res2 = list_connected_partitions(4);
    
//    for( const auto& it : res2 )
//    {
//        for( const auto& inner : it )
//        {
//            for( const auto& inner2 : inner )
//            {
//                std::cout << inner2 << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
    
    return 0;
}
