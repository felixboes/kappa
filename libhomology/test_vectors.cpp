#include <iostream>

#include "homology.hpp"

int main()
{
    VectorQ v(3);
    v(0) = 1;
    v(1) = 2;
    v(2) = -3;
    
    std::cout << v << std::endl;
    
    VectorBool w(3);
    w.add_entry(0);
    w.add_entry(2);
    
    std::cout << w << std::endl;
    
    v.resize(10);
    std::cout << v << std::endl;
    return 0;
}
