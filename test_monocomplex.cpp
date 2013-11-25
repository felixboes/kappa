#include <iostream>

#include <homology.hpp>
#include "monocomplex.hpp"


int main()
{
    typedef MonoComplex<Q> MonoComplexQ;
    MonoComplexQ mc(1,2);
    mc.show_basis(1);
    mc.show_basis(2);
    mc.show_basis(3);
    mc.show_basis(4);
    mc.show_basis(5);
    
    mc.gen_differential(4);
    
    return 0;
}
