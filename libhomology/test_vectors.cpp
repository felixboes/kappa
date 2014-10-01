#include <iostream>

#include "homology.hpp"

int main()
{
    VectorQ v(3);
    v(0) = 1;
    v(1) = 2;
    v(2) = -3;
    
    VectorQ v_2(3);
    v_2(0) = 3;
    v_2(1) = Q(7) / Q(4);
    v_2(2) = -2;
    
    v += v_2;
    
    std::cout << v << std::endl;
    
    VectorBool w(3);
    w.add_entry(0);
    w.add_entry(2);
    
    VectorBool w_2(3);
    w_2.add_entry(0);
    w_2.add_entry(1);
    
    w += w_2;
    std::cout << w << std::endl;
    
    v.resize(10);
    std::cout << v << std::endl;
    
    w.resize(10);
    std::cout << w << std::endl;
    
    return 0;
}
