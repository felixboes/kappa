
#include <iostream>
#include <tuple>
#include <type_traits>

#include "kappa.hpp"

template < size_t n, typename... T >
typename std::enable_if< ( n >= sizeof...(T) ) >::type print_tuple( std::ostream&   , const std::tuple< T... >& )
{}

template < size_t n, typename... T >
typename std::enable_if< ( n  < sizeof...(T) ) >::type print_tuple( std::ostream& os, const std::tuple< T... >& tup )
{
    if ( n != 0 )
    {
        os << ", ";
    }
    os << std::get<n>( tup );
    print_tuple< n+1 >( os, tup );
}

template < typename... T >
std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup)
{
    os << "[";
    print_tuple< 0 >( os, tup );
    return os << "]";
}

int main( int , char**  )
{
    OperationTesterQ Opt;
    OperationTesterQ::MonoIndex idx_v( true, 1, 2, 8 );
    OperationTesterQ::MonoIndex idx_w( true, 2, 0, 7 );
    
    std::cout << idx_v << std::endl;
    std::cout << idx_w << std::endl;
    std::cout << OperationTesterQ::product(idx_v, idx_w) << std::endl;
    
    Opt.load_basis(idx_v);
    Opt.load_basis(true, 0, 2, 3, true);
    std::cout << Opt.dim( true, 0, 2, 3 ) << std::endl;
    std::cout << Opt.dim( true, 1, 2, 7 ) << std::endl;

    Opt.forget_basis(true, 0, 2, 3);
    std::cout << Opt.dim( true, 0, 2, 3 ) << std::endl;
    
    return 0;
}
