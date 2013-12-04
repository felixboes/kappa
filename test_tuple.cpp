#include "factorial.hpp"
#include "tuple.hpp"

#include <iostream>

bool compare_boundary_computation (Tuple t, uint32_t i )
{
    Tuple boundary_1 = t.d_hor(i);
    Tuple boundary_2 = t.d_hor_naive(i);
    std::cout << boundary_1 << " " << boundary_2 << ( boundary_1 == boundary_2? " test passed" : " test failed" ) << std::endl;
    
    return boundary_1 == boundary_2;
}

bool compare_boundary_computations (Tuple t )
{
    bool ret = false;
    for( uint32_t i = 0; i <= t.p ; ++i)
    {
        ret = compare_boundary_computation(t,i);
    }
    return ret;
}

int main()
{
    std::cout.setf(std::ios::unitbuf);
    
    Tuple t(8, 5);
    
    t[1] = Transposition(5,1);
    t[2] = Transposition(7,4);
    t[3] = Transposition(8,6);
    t[4] = Transposition(8,3);
    t[5] = Transposition(8,2);
    
    std::cout << "t: " << t << std::endl;
    std::cout << "n_c n_p:" << t.permutation_type().num_cycles << " " << t.permutation_type().num_punctures << std::endl;
    std::cout << "m" << t.monotone() << std::endl;
    
    std::cout << "t: " << t << std::endl;
    
    std::vector<uint8_t> m = {3,2,4,3};
    for( auto& it : m )
    {
        std::cout << (uint32_t)it << " " << (  t.f(it)?"not degenerate":"    degenerate" ) << " " << t << std::endl;
    }
    
    t[1] = Transposition(5,1);
    t[2] = Transposition(7,4);
    t[3] = Transposition(8,6);
    t[4] = Transposition(8,3);
    t[5] = Transposition(8,2);
    
    compare_boundary_computations(t);
    t.phi(1,1);
    compare_boundary_computations(t);
    t.phi(2,2);
    compare_boundary_computations(t);
    t.phi(3,3);
    compare_boundary_computations(t);
    t.phi(4,2);
    compare_boundary_computations(t);
    t.phi(5,3);
    compare_boundary_computations(t);
    
    HashTuple h;
    std::cout << "Hash of " << t << ":" << h(t) << std::endl;
    
    std::cout << t << std::endl;
}
