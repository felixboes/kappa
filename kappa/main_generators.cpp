// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#include "kappa.hpp"

//template< class CoefficientT >
//void cohomology_generators( const uint32_t g, const uint32_t m, const int32_t p )
//{
//    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
//    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    MatrixField< CoefficientT > vanishing_test = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    MonoBasis basis = load_parallel_mono_basis( g, m, p );

//    auto base = compute_base_of_kernel< MatrixField<CoefficientT>, VectorField<CoefficientT> >( kernel );
//    for( const auto& v : base )
//    {
//        auto c = v.homology_class( kernel, image );
//        if( c.is_zero() == false )
//        {
//            std::cout << "Chain = ";

//            for( const auto& it : basis.basis )
//            {
//                if( v.at(it.id) != CoefficientT(0) )
//                {
//                    std::cout << std::setw(4) << v.at(it.id) << " * " << it << std::endl;
//                }
//            }

//            std::cout << "The number of affiliated cells is " << v.number_non_vanishing_entries() << std::endl
//                      << "This cochain is " << ( matrix_vector_product_vanishes(vanishing_test, v) == true ? "indeed " : "NOT " ) << "a cocycle." << std::endl
//                      << "Its class it " << c << std::endl
//                      << std::endl;
//        }
//    }
//}

//template< class CoefficientT >
//void cohomology_generators_tex( const uint32_t g, const uint32_t m, const int32_t p )
//{
//    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p-1) + "_base_changes", false );
//    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    MatrixField< CoefficientT > vanishing_test = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_parallel_differentials<CoefficientT>(g,m) + std::to_string(p) + "_triangular", false );
//    MonoBasis basis = load_parallel_mono_basis( g, m, p );


//    std::cout << tex_preamble();
//    auto base = compute_base_of_kernel< MatrixField<CoefficientT>, VectorField<CoefficientT> >( kernel );
//    DoubleComplex< ChainComplex< CoefficientT, MatrixField<CoefficientT>, DiagonalizerField<MatrixField<CoefficientT>>, HomologyField > > dc(g, m, SignConvention::all_signs, 1, 1);
//    for( const auto& v : base )
//    {
//        auto c = v.homology_class( kernel, image );
//        if( c.is_zero() == false )
//        {
//            std::cout << "The following cochain represents the class " << c << " $\\in H_" << (4*g + 2*m - p) << "(\\mathfrak M_{" << g << ",1}^{" << m << "};\\mathbb K)$\\\\" << std::endl;
//            dc.compute_proj_E(p+1);

//            for( const auto& it : basis.basis )
//            {
//                if( v.at(it.id) != CoefficientT(0) )
//                {
//                    dc.proj_E_ast_tex(v.at(it.id), it);
//                }
//            }

//            std::cout << "\\clearpage" << std::endl;
//        }
//    }
//    std::cout << tex_end();
//}

//template< class CoefficientT >
//void print_cohomology_generators_tex( const uint32_t g=2, const uint32_t m=0, const int32_t p=3 )
//{
//    cohomology_generators_tex<CoefficientT>(g, m, 4*g+2*m-p);
//}

template< class CoefficientT >
VectorField< CoefficientT > cohomology_class( const uint32_t g, const uint32_t m, const int32_t p, const bool radial_model_used, const VectorField< CoefficientT >& v )
{         
    MatrixField< CoefficientT > image  = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_differentials<CoefficientT>(g,m,radial_model_used) + std::to_string(p-1) + "_base_changes", false );
    MatrixField< CoefficientT > kernel = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_differentials<CoefficientT>(g,m,radial_model_used) + std::to_string(p) + "_triangular", false );
//    std::cout << "Image:" << std::endl << image << std::endl << std::endl;
//    std::cout << "Kernel:" << std::endl << kernel << std::endl << std::endl;
//    std::cout << "Diagonal:" << std::endl << diagonal << std::endl << std::endl;
    const MatrixField< CoefficientT > M = load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_differentials<CoefficientT>(g,m,radial_model_used) + std::to_string(p) + "_triangular", false );
    const bool v_is_a_cocycle = matrix_vector_product_vanishes(M, (VectorField<CoefficientT>)v);
    if( v_is_a_cocycle )
    {
        std::cout << "This cochain is indeed a cocycle." << std::endl;
    }
    else
    {
        std::cout << "This cochain is NOT a cocycle. See: ";
        const auto w = matrix_vector_product(load_from_file_bz2< MatrixField< CoefficientT > >( filename_prefix_differentials<CoefficientT>(g,m,radial_model_used) + std::to_string(p) + "_triangular", false ), v );
        std::cout << w << std::endl;
    }
    return v.homology_class( kernel, image );
}

template< class CoefficientT >
void test_( const MonoCochainField< CoefficientT >& cochain )
{
    std::cout << "Name                  = " << cochain.get_name() << std::endl;
    std::cout << "(g,m,p)               = (" << cochain.get_g() << "," << cochain.get_m() << "," << 4*cochain.get_g() + 2*cochain.get_m() - cochain.get_p() << ")" << std::endl;
    std::cout << "Cohomology class      = " << cohomology_class( cochain.get_g(), cochain.get_m(), cochain.get_p(), cochain.get_radial(), cochain ) << std::endl;
    std::cout << std::endl;
}

template< class CoefficientT >
void test_( const std::string& name, const uint32_t g, const uint32_t m, const uint32_t homological_p, const Tuple& cell )
{
    MonoCochainField< CoefficientT > cochain(g, m, 4*g+2*m-homological_p );
    cochain.add_kappa_dual( CoefficientT(1),cell);
    
    test_< CoefficientT >(name, g, m, homological_p, cochain);
}

void verify_known_generators()
{
    std::cout << "Mod 2 computations." << std::endl;
    std::cout << "--------------------------------" << std::endl;
    Zm::set_modulus(2);

    // Classes for g = 0:
    // m = 1
    test_( create_cochain<Zm>( Generator::a ) );
    // m = 2
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) );
    test_( create_cochain<Zm>( Generator::b ) );
    // m = 4
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) );
    test_( create_cochain<Zm>( Generator::Qb ) );
    // m = 6
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Qb ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::Qb ) );

    // Classes for g = 1:
    // m = 0:
    test_( create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::d ) );
    // m = 1:
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::Eb ) );
    // m = 2:
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::f ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
    // m = 3:
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::Eb ) );
    // m = 4:
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::e ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::Eb ) );
    test_( create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::a ) * create_cochain<Zm>( Generator::b ) * create_cochain<Zm>( Generator::Eb ) );
     
    // Classes for g = 2:
    // m = 0:
    test_( create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::c ) );
    test_( create_cochain<Zm>( Generator::c ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::Q_alpha_inv_c ) );
    test_( create_cochain<Zm>( Generator::Q_beta_c ) );
    test_( create_cochain<Zm>( Generator::Q_gamma_c ) );
    test_( create_cochain<Zm>( R_alpha_inv_beta_c_c ) );
    test_( create_cochain<Zm>( Generator::d ) * create_cochain<Zm>( Generator::d ) );
    test_( create_cochain<Zm>( Generator::Te ) );
    test_( create_cochain<Zm>( Generator::Q_alpha_inv_d ) );
    test_( create_cochain<Zm>( Generator::Q_beta_d ) );
    test_( create_cochain<Zm>( Generator::Q_gamma_d ) );
    test_( create_cochain<Zm>( R_alpha_inv_beta_c_d ) );
    test_( create_cochain<Zm>( R_alpha_inv_gamma_inv_c_d ) );
    test_( create_cochain<Zm>( R_alpha_inv_gamma_c_d ) );
    test_( create_cochain<Zm>( R_alpha_inv_beta_inv_c_d ) );
    test_( create_cochain<Zm>( R_alpha_inv_alpha_inv_c_d ) );
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

int main( int argc, char** argv )
{
    std::cout.setf(std::ios::unitbuf);
    
    std::cout << kappa_version( argc, argv ) << std::endl;
    
    if( argc == 1 )
    {
        verify_known_generators();
    }
    else
    {

//      std::cout << "Rational computations." << std::endl;
//      std::cout << "--------------------------------" << std::endl;
//      test_( create_cochain<Q>( Generator::R_alpha_inv_alpha_inv_c_d ) );

//      std::cout << "Mod 2 computations." << std::endl;
//      std::cout << "--------------------------------" << std::endl;
//      Zm::set_modulus(2);
//      test_( create_cochain<Zm>( Generator::radial_a ) );
//      test_( create_cochain<Zm>( Generator::radial_b ) );
//      test_( create_cochain<Zm>( Generator::radial_Srad_Te) );
//      test_( create_cochain<Zm>( Generator::T_1_f ) );
//      test_( create_cochain<Zm>( Generator::T_2_f ) );

//      std::cout << "Mod 5 computations." << std::endl;
//      std::cout << "--------------------------------" << std::endl;
//      Zm::set_modulus(5);
//      test_( create_cochain<Zm>( Generator::Q_beta_c ) );
//      test_( create_cochain<Zm>( Generator::Q_gamma_c ) );
//      test_( create_cochain<Zm>( Generator::R_alpha_inv_alpha_inv_c_d ) );
    }

    return 0;
}
