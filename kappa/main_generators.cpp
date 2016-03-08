#include "kappa.hpp"

template< class CoefficientT >
VectorField< CoefficientT > cohomology_class( const uint32_t g, const uint32_t m, const int32_t p, const VectorField< CoefficientT > v )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    std::cout << "Image:" << std::endl << image << std::endl << std::endl;
//    std::cout << "Kernel:" << std::endl << kernel << std::endl << std::endl;
//    std::cout << "Diagonal:" << std::endl << diagonal << std::endl << std::endl;
    const bool v_is_a_cocycle = matrix_vector_product_vanishes(load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false ), v);
    if( v_is_a_cocycle )
    {
        std::cout << "This cochain is indeed a cocycle." << std::endl;
    }
    else
    {
        std::cout << "This cochain is NOT a cocycle. See: ";
        const auto w = matrix_vector_product(load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false ), v );
        std::cout << w << std::endl;
    }
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
                      << "This cochain is " << ( matrix_vector_product_vanishes(vanishing_test, v) == true ? "indeed " : "NOT " ) << "a cocycle." << std::endl
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
    DoubleComplex< ChainComplex< CoefficientT, MatrixField<CoefficientT>, DiagonalizerField<MatrixField<CoefficientT>>, HomologyField > > dc(g, m, SignConvention::all_signs, 1, 1);
    for( const auto& v : base )
    {
        auto c = v.homology_class( kernel, image );
        if( c.is_zero() == false )
        {
            std::cout << "The following cochain represents the class " << c << " $\\in H_" << (4*g + 2*m - p) << "(\\mathfrak M_{" << g << ",1}^{" << m << "};\\mathbb K)$\\\\" << std::endl;
            dc.compute_proj_E(p+1);
            
            for( const auto& it : basis.basis )
            {
                if( v.at(it.id) != CoefficientT(0) )
                {
                    dc.proj_E_ast_tex(v.at(it.id), it);
                }
            }
            
            std::cout << "\\clearpage" << std::endl;
        }
    }
    std::cout << tex_end();
}

template< class CoefficientT >
void test_( const MonoCochainField< CoefficientT >& cochain )
{
    std::cout << "Name                  = " << cochain.get_name() << std::endl;
    std::cout << "(g,m,p)               = (" << cochain.get_g() << "," << cochain.get_m() << "," << 4*cochain.get_g() + 2*cochain.get_m() - cochain.get_p() << ")" << std::endl;
    std::cout << "Cohomology class      = " << cohomology_class( cochain.get_g(), cochain.get_m(), cochain.get_p(), cochain ) << std::endl;
    std::cout << std::endl;
}

template< class CoefficientT >
void test_( const std::string& name, const uint32_t g, const uint32_t m, const uint32_t homological_p, const Tuple& cell )
{
    MonoCochainField< CoefficientT > cochain(g, m, 4*g+2*m-homological_p );
    cochain.add_kappa_dual( CoefficientT(1),cell);
    
    test_< CoefficientT >(name, g, m, homological_p, cochain);
}

template< class CoefficientT >
MonoCochainField< CoefficientT > create_z_1()
{
    MonoCochainField< CoefficientT > cochain(2, 0, 5);
    cochain.set_name("z");
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 2, 3, 1, 2, 1 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 2, 3, 1, 2, 1 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 4, 2, 3, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 4, 2, 3, 2 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 2, 4, 1, 4, 2, 3, 2 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 4, 2, 3, 2 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 2, 1, 4, 2, 3, 2 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 2, 1, 4, 2, 3, 2 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 3, 1, 4, 2, 3, 2 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 4, 2, 3, 2, 5, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 4, 1, 3, 1, 5, 2, 3, 2 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
    
    return cochain;
}

template< class CoefficientT >
void test_z_1()
{
    test_( create_z_1<CoefficientT>() );
}

template< class CoefficientT >
MonoCochainField< CoefficientT > create_z_2()
{
    MonoCochainField< CoefficientT > cochain(2, 0, 5);
    cochain.set_name("z");
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 3, 3, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 1, 3, 1, 2, 1 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 1, 3, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 1, 4, 3, 4, 1, 2, 1 ) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 4, 3, 4, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 5, 3, 3, 1, 4, 1, 2, 1 ) );
    
    cochain.add_kappa_dual( CoefficientT(1), create_cell( 4, 4, 1, 2, 1, 5, 3, 4, 3 ) );
    
    return cochain;
}

template< class CoefficientT >
void test_z_2()
{
    test_( create_z_2<CoefficientT>() );
}

template< class CoefficientT >
void test_tilde_d()
{
    MonoCochainField< CoefficientT > cochain(2, 0, 5);
    cochain.set_name("tilde d");
    cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 4, 5, 4, 3, 1, 2, 1) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 3, 4, 3, 2, 1, 2, 1) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 4, 3, 5, 4, 2, 1, 2, 1) );
    cochain.add_kappa_dual( CoefficientT(1), create_cell(4, 5, 4, 5, 4, 2, 1, 3, 2) );
    test_( cochain );
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

void verify_known_generators()
{
    std::cout << "Mod 2 computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    Zm::set_modulus(2);

//     // Classes for g = 0:
//     // m = 1
//     test_( create_cochain<Zm>( Generator::a ) );
//     // m = 2
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) );
//     test_( create_cochain<Zm>( Generator::b ) );
//     // m = 4
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) );
//     test_( create_cochain<Zm>( Generator::Qb ) );
//     // m = 6
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Qb ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::Qb ) );
//     
//     // Classes for g = 1:
//     // m = 0:
//     test_( create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::d ) );
//     // m = 1:
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::e ) );
//     test_( create_cochain<Zm>( Generator::Eb ) );
//     // m = 2:
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::e ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
//     // m = 3:
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) ); 
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::e ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::e ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::Eb ) );
//     // m = 4:
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) ); 
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::e ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::e ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
//     test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
//     test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::Eb ) );
//     
    // Classes for g = 2:
    // m = 0:
    test_( create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::d ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::Qd ) );
    test_( create_cochain<Zm>( Generator::Te ) );
    test_( create_cochain<Zm>( Generator::TEb ) );
    // m = 1:
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Qd ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Te ) );
    test_( create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::Eb ) );
    test_( create_cochain<Zm>( Generator::d ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::TEb ) );
    test_( create_cochain<Zm>( Generator::d ) * create_cochain<Zm>( Generator::Eb ) );
}

template< class CoefficientT >
void print_cohomology_generators_tex( const uint32_t g=2, const uint32_t m=0, const int32_t p=3 )
{
    cohomology_generators_tex<CoefficientT>(g, m, 4*g+2*m-p);
}

int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << kappa_version( argc, argv ) << std::endl;
    
    if( argc > 1 )
    {
        if( atoi( argv[1] ) == 1 )
        {
            verify_known_generators();
            return 0;
        }
        else if( atoi(argv[1]) == 2 )
        {
            if( argc > 4 )
            {
                print_cohomology_generators_tex<Q>( atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) );
            }
            else
            {
                print_cohomology_generators_tex<Q>();
            }
            return 0;
        }
    }
    
    std::cout << "Rational computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;

    test_( create_cochain<Q>( Generator::R_alpha_inv_alpha_inv_c_d ) );
    test_( create_cochain<Q>( Generator::u) );

    std::cout << "Mod 2 computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    Zm::set_modulus(2);

    // Moreover:
//    test_( create_cochain<Zm>( Generator::Qc ) );
//    test_( create_cochain<Zm>( Generator::Q_alpha_c ) );
//    test_( create_cochain<Zm>( Generator::Q_beta_c ) );
//    test_( create_cochain<Zm>( Generator::Q_gamma_c ) );
    
//    test_( create_cochain<Zm>( Generator::Te ) );
//    test_( create_cochain<Zm>( Generator::Q_alpha_d ) );
//    test_( create_cochain<Zm>( Generator::Q_beta_d ) );
//    test_( create_cochain<Zm>( Generator::Q_gamma_d ) );
    
//     test_( create_cochain<Zm>( Generator::R_alpha_inv_alpha_inv_c_d ) );
    
    // Later:
//    test_z_1<Zm>();
//    test_z_2<Zm>();
//    test_tilde_d<Zm>();
    
    
//     test_( create_cochain<Zm>( Generator::R_a_e) );
    
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
    test_( create_cochain<Zm>( Generator::u) );
    
    std::cout << "Mod 5 computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    Zm::set_modulus(5);
    
//    test_( create_cochain<Zm>( Generator::Q_alpha_c ) );
//    test_( create_cochain<Zm>( Generator::Q_beta_c ) );
//    test_( create_cochain<Zm>( Generator::Q_gamma_c ) );
//     test_( create_cochain<Zm>( Generator::R_alpha_inv_alpha_inv_c_d ) );
    
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
