#include "factorial.hpp"
#include "tuple.hpp"

#include <iostream>

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
    std::cout << "m" << t.monoton() << std::endl;
    
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
    
    
    t.phi(1,1);
    t.phi(2,2);
    t.phi(3,3);
    t.phi(4,2);
    t.phi(5,3);
    
    for( uint8_t i = 0; i <= 8; ++i)
    {
        std::cout << t.d_hor(i) << " " << t.d_hor_test(i) << (t.d_hor(i) == t.d_hor_test(i)? " test passed" : " test failed" ) << std::endl;
        
    }
        
    std::cout << t << std::endl;
}
