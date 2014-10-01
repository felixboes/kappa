#include <iostream>

#include "homology.hpp"

int main()
{
    VectorQ v(3);
    v(0) = 1;
    v(1) = 0;
    v(2) = 0;
    
    // 2,  1,  1
    // 1, -1,  2
    // 0,  1, -1
    MatrixQ m(3,3);
    m(0,0) = 2;
    m(0,1) = 1;
    m(0,2) = 1;
    m(1,0) = 1;
    m(1,1) = -1;
    m(1,2) = 2;
    m(2,0) = 0;
    m(2,1) = 1;
    m(2,2) = -1;
    
    std::cout << m << std::endl;
    DiagonalizerField<MatrixQ> diago;
    diago(m);
    std::cout << m << std::endl;
    
    VectorQ v_prime(3);
    v_prime(0) = 0;
    v_prime(1) = 1;
    v_prime(2) = 1;
    apply_base_changes( m, v_prime );
    std::cout << v_prime << std::endl;
    
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
