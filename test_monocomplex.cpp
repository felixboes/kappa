#include <iostream>

#include <homology.hpp>
#include "monocomplex.hpp"


int main()
{
    typedef MonoComplex<ChainComplexQ> MonoComplexQ;
    MonoComplexQ mc(1,2);
//    mc.show_basis(1);
//    mc.show_basis(2);
//    mc.show_basis(3);
//    mc.show_basis(4);
//    mc.show_basis(5);
    
//    mc.gen_differential(4);
    
    mc.gen_differential(4);
    mc.gen_differential_naive(4);
    
    mc.show_differential(4);
    mc.show_differential_naive(4);
    
    return 0;
}
