#include <iostream>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures (rational:0 | Z_r: r > 0)" << std::endl;
}

int main(int argc, char** argv)
{
    SessionConfig conf(argc, argv);
    if( conf.valid == false )
    {
        print_usage(argc, argv);
        return 1;
    }
    
    // We may start with the computations.
    if(conf.rational == true)
    {
        MonoComplexQ monocomplex( conf.genus, conf.num_punctures );
        monocomplex.gen_differentials();
        
        for( auto& it : monocomplex.basis_complex )
        {
            auto& p = it.first;
            monocomplex.show_basis(p);
            monocomplex.show_differential(p);
        }
        
        std::cout << monocomplex.matrix_complex.homology() << std::endl;
        
        return 0;
    }
    else
    {
        MonoComplexZm monocomplex( conf.genus, conf.num_punctures );
        monocomplex.gen_differentials();
        
//        for( auto& it : monocomplex.basis_complex )
//        {
//            auto& p = it.first;
//            monocomplex.show_basis(p);
//            monocomplex.show_differential(p);
//        }
        
        std::cout << monocomplex.matrix_complex.homology() << std::endl;
        
        return 0;
    }
}
