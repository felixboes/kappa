
#include <iostream>
#include <tuple>
#include <type_traits>

#include "kappa.hpp"

int main( int , char**  )
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
    Opt.load_basis(true, 1, 2, 3);
    Opt.load_basis(true, 1, 2, 5);
    Opt.load_basis(true, 1, 2, 6);
    std::cout << Opt.dim( true, 0, 2, 3 ) << std::endl;
    std::cout << Opt.dim( true, 1, 2, 7 ) << std::endl;
    if( Opt.dim(true, 0, 2, 3) > 0 )
    {
        std::cout << Opt.basis.at( MonoIndex(true, 0, 2, 3) ) << std::endl;
    }
    
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
    
    Opt.forget_base_changes(true, 2, 2, 5);
    Opt.load_triangular(true, 2, 2, 5);
    Opt.print_cache_status();
    
    return 0;
}
