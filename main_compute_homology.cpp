#include <iostream>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures (rational:0 | Z_r: r > 0)" << std::endl;
}

int main(int argc, char** argv)
{
    // Parse configuration from command line arguments.
    SessionConfig conf(argc, argv);
    if( conf.valid == false )
    {
        print_usage(argc, argv);
        return 1;
    }
    
    if ( conf.setup_configuration() == false )
    {
        std::cout << "The configuration could not been setup." << std::endl;
        return 2;
    }
    
    // We may start with the computations.
    if(conf.rational == true)
    {
        // Compute all bases.
        MonoComplexQ monocomplex( conf.genus, conf.num_punctures );
        // Compute all differentials.
        monocomplex.gen_differentials();
        
//        // Print all bases and differentials to screen.
//        for( auto& it : monocomplex.basis_complex )
//        {
//            auto& p = it.first;
//            monocomplex.show_basis(p);
//            monocomplex.show_differential(p);
//        }
        
        // Compute homology and print to screen.
        std::cout << monocomplex.matrix_complex.homology() << std::endl;
        
        return 0;
    }
    else
    {
        // Compute all bases.
        MonoComplexZm monocomplex( conf.genus, conf.num_punctures );
        // Compute all differentials.
        monocomplex.gen_differentials();
        
//        // Print all bases and differentials to screen.
//        for( auto& it : monocomplex.basis_complex )
//        {
//            auto& p = it.first;
//            monocomplex.show_basis(p);
//            monocomplex.show_differential(p);
//        }
        
        // Compute homology and print to screen.
        std::cout << monocomplex.matrix_complex.homology() << std::endl;
        
        return 0;
    }
}
