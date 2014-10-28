
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

void test_matrix_vector_stuff()
{
    OperationTesterQ Opt("q");
    OperationTesterQ::MonoIndex idx(true, 0, 2, 4);
    Opt.load_basis(idx);
    std::cout << Opt.basis[idx] << std::endl;
    
    Opt.load_diagonal(idx);
    std::cout << Opt.diagonal[idx] << std::endl;
    
    Opt.load_base_changes(idx);
    std::cout << Opt.base_changes[idx] << std::endl;
    
    Opt.load_triangular(idx);
    std::cout << Opt.triangular[idx] << std::endl;
    
    VectorQ v(2);
    v(0) = 1;
    
    VectorQ w(2);
    w(1) = 1;
    
    std::cout << v << std::endl;
    std::cout << w << std::endl << std::endl;
    
    apply_base_changes( Opt.base_changes[idx], v );
    apply_base_changes( Opt.base_changes[idx], w );
    
    std::cout << v << std::endl;
    std::cout << w << std::endl << std::endl;
    
    std::cout << Opt.triangular[idx] << std::endl;
    std::cout << Opt.diagonal[idx] << std::endl;
    
    VectorQ c(2);
    c(0) = 4;
    c(1) = 4;
    std::cout << "The vector " << c << " is " << (Opt.vector_is_cycle( idx, c ) == true ? "indeed " : "not " ) << "a cycle." << std::endl;
    std::cout << matrix_vector_product( Opt.triangular[idx], c ) << std::endl;
    
    Opt.forget_triangular(idx);
    
    VectorQ d(2);
    d(0) = 1;
    d(1) = 2;
    std::cout << "The vector " << c << " is " << (Opt.vector_is_cycle( idx, d ) == true ? "indeed " : "not " ) << "a cycle." << std::endl;
    std::cout << matrix_vector_product( Opt.triangular[idx], d ) << std::endl;
    
    Opt.vector_print_homology_class( idx, c );
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
            test_matrix_vector_stuff();
        }
        else
        {
            goto remaining_stuff_to_do;
        }
        return 0;
    }
    remaining_stuff_to_do:
    
    test_matrix_vector_stuff();
    
    return 0;
}
