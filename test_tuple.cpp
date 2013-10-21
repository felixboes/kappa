#include "factorial.hpp"
#include "tuple.hpp"

#include <iostream>

int main()
{
    std::cout.setf(std::ios::unitbuf);
    
    Tuple t(5);
    t.p = 8;
    
    t[0] = Transposition(5,1);
    t[1] = Transposition(7,4);
    t[2] = Transposition(8,6);
    t[3] = Transposition(8,3);
    t[4] = Transposition(8,2);
    
    std::cout << "t: " << t << std::endl;
    std::cout << "n_c n_p:" << t.permutation_type().num_cycles << " " << t.permutation_type().num_punctures << std::endl;
    std::cout << "m" << t.monoton() << std::endl;
    
    std::vector<uint8_t> m = {1, 2, 3 ,2,3 ,3,4};
    for( auto& it : m )
    {
        std::cout << (uint32_t)it << " " << t.f(it) << " " << t << std::endl;
    }
}
