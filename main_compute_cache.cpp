#include <iostream>
#include <string>

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
        MonoComplexZStorageOnly monocomplex( conf.genus, conf.num_punctures );
        
        // Save bases to file
        std::string prefix_basis("./cache/bases/");
        prefix_basis += std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_";
        
        for( auto& it : monocomplex.basis_complex )
        {
            // Store the p-th basis.
            auto& p = it.first;
            save_to_file_bz2<MonoBasis>(it.second, prefix_basis + std::to_string(p));
        }
        
        // Compute all differentials.
        std::string prefix_differentials("./cache/differentials/");
        prefix_differentials += std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_";
        for( auto& it : monocomplex.basis_complex )
        {
            auto& p = it.first;
            monocomplex.gen_differential( p );
            save_to_file_bz2<MatrixZDontDiagonalize>( monocomplex.matrix_complex[p], prefix_differentials + std::to_string(p) );
            monocomplex.delete_differential(p);
        }
        
        return 0;
    }
}
