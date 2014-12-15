
#include <iostream>
#include <tuple>
#include <type_traits>

#include "kappa.hpp"

void test_member_methods()
{
    OperationTesterQ Opt( "q" );
    typedef OperationTesterQ::MonoIndex MonoIndex;
    OperationTesterQ::MonoIndex idx_v( true, 1, 2, 8 );
    OperationTesterQ::MonoIndex idx_w( true, 2, 0, 7 );
    
    std::cout << idx_v << std::endl;
    std::cout << idx_w << std::endl;
    std::cout << OperationTesterQ::product(idx_v, idx_w) << std::endl;
    
    Opt.load_basis(idx_v);
    Opt.load_basis(true, 0, 2, 3, true);
    Opt.load_basis(true, 0, 1, 2);
    Opt.load_basis(true, 0, 2, 2);
    Opt.load_basis(true, 0, 2, 3);
    Opt.load_basis(true, 0, 2, 4);
    Opt.load_basis(true, 1, 2, 3);
    Opt.load_basis(true, 1, 2, 5);
    Opt.load_basis(true, 1, 2, 6);
    std::cout << Opt.dim( true, 0, 2, 3 ) << std::endl;
    std::cout << Opt.dim( true, 1, 2, 7 ) << std::endl;
    if( Opt.dim(true, 0, 2, 3) > 0 )
    {
        std::cout << Opt.basis.at( MonoIndex(true, 0, 2, 3) ) << std::endl;
    }
    std::cout << Opt.basis.at( MonoIndex(true, 0, 1, 2) ) << std::endl;
    std::cout << Opt.basis.at( MonoIndex(true, 0, 2, 2) ) << std::endl;
    std::cout << Opt.basis.at( MonoIndex(true, 0, 2, 3) ) << std::endl;
    std::cout << Opt.basis.at( MonoIndex(true, 0, 2, 4) ) << std::endl;
    
    Opt.forget_basis(true, 0, 2, 3);
    std::cout << Opt.dim( true, 0, 2, 3 ) << std::endl;
    
    Opt.load_base_changes(true, 0,1,2);
    Opt.load_base_changes(true, 1,2,4);
    Opt.load_base_changes(true, 2,2,5);
    Opt.forget_base_changes(true, 0,1,2);
    Opt.forget_base_changes(true, 0,2,4);
    
    Opt.load_triangular(true, 0,1,2);
    Opt.load_triangular(true, 0,2,4);
    Opt.forget_triangular(true, 0,2,4);
    
    Opt.load_diagonal(true, 0,1,2);
    Opt.load_diagonal(true, 0,2,4);
    Opt.forget_diagonal(true, 0,1,2);
    
    Opt.print_cache_status();
    
    std::cout << Opt.base_changes[MonoIndex(true,2,2,5)].diagonal << std::endl;
    
    Opt.forget_base_changes(true, 2, 2, 5);
    Opt.load_triangular(true, 2, 2, 5);
    Opt.print_cache_status();
    
    Opt.load_diagonal(true, 2, 2, 5);
    std::cout << Opt.diagonal[MonoIndex(true,2,2,5)] << std::endl;
}

void test_products()
{
    OperationTesterQ Opt("q");
    OperationTesterQ::MonoIndex idx(true, 0, 2, 4);
    OperationTesterQ::MonoIndex idx_prod(true, 0, 4, 8);
    Opt.load_basis(idx);
    Opt.load_basis(idx_prod);
    
    std::cout << Opt.basis[idx] << std::endl;
    std::cout << Opt.basis[idx_prod] << std::endl;
    
    Tuple t_1 = *( Opt.basis[idx].basis.begin() );
    Tuple t_2 = *( std::next( Opt.basis[idx].basis.begin() ) );
    Tuple t_prod = t_1*t_2;
    t_prod.id = Opt.basis[idx_prod].id_of( t_prod );
    std::cout << t_2 << " * " << t_1 << " = " << t_prod << std::endl;
    
    VectorQ v_1 (2);
    v_1(0) = 2;
    v_1(1) = 5;
    
    VectorQ v_2 (2);
    v_2(0) = 0;
    v_2(1) = 7;
       
    std::cout << v_1 << " * " << v_2 << " = " << Opt.product(idx, v_1, idx, v_2) << std::endl;
}

void test_kappa_dual()
{
    OperationTesterQ Opt("q");
    OperationTesterQ::MonoIndex idx(true, 0, 1, 2);
    OperationTesterQ::MonoIndex idx_res(true, 0, 2, 4);
    Opt.load_basis( idx );
    Opt.load_basis( idx_res );
    
    Tuple t(2);
    t[1] = Transposition( 4, 3 );
    t[2] = Transposition( 2, 1 );
    std::cout << t << std::endl;
    
    VectorQ v(2);
    std::vector< size_t > s;
    s.push_back(1);
   
    Q c(1);
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t[1] = Transposition( 4, 2 );
    t[2] = Transposition( 3, 1 );
    std::cout << t << std::endl;
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t[1] = Transposition( 4, 1 );
    t[2] = Transposition( 3, 2 );
    std::cout << t << std::endl;
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t[1] = Transposition( 3, 2 );
    t[2] = Transposition( 3, 1 );
    std::cout << t << std::endl;
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t[1] = Transposition( 3, 2 );
    t[2] = Transposition( 2, 1 );
    std::cout << t << std::endl;
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t[1] = Transposition( 3, 1 );
    t[2] = Transposition( 2, 1 );
    std::cout << t << std::endl;
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t[1] = Transposition( 2, 1 );
    t[2] = Transposition( 3, 1 );
    std::cout << t << std::endl;
    Opt.compute_and_add_kappa_dual_rec(c, t, Opt.basis.at(idx_res), v, s, 0);
    
    t = Tuple(3);
    t[1] = Transposition( 5, 1 );
    t[2] = Transposition( 3, 2 );
    t[3] = Transposition( 4, 3 );
    Opt.compute_and_add_kappa_dual(c, t, Opt.basis.at(idx_res), v);
}

int main( int argc , char** argv )
{
    if( argc > 1 )
    {
        if( atoi(argv[1]) == 1)
        {
            test_member_methods();
        }
        else if ( atoi(argv[1]) == 2 )
        {
            test_products();
        }
        else if( atoi(argv[1]) == 3 )
        {
            test_kappa_dual();
        }
        else
        {
            goto remaining_stuff_to_do;
        }
        return 0;
    }
    remaining_stuff_to_do:
    
    test_products();
    
    return 0;
}
