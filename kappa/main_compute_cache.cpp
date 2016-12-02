#include <iostream>
#include <string>

#include <libhomology/homology.hpp>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg" << std::endl;
}

void initialize_zero_matrix( MatrixZDontDiagonalize & matrix, size_t num_rows, size_t num_cols)
{
    matrix.resize(num_rows, num_cols);
    for (size_t i = 0; i < num_rows; ++i)
    {
        for (size_t j = 0; j < num_cols; ++j)
        {
            matrix(i, j) = 0;
        }
    }
}

int main(int argc, char** argv)
{
    // Parse configuration from command line arguments.
    SessionConfig conf(argc, argv);
    if( conf.option_set("help") )
    {
        print_usage(argc, argv);
        std::cout << conf.desc << std::endl;
        return 0;
    }
    if( ! ( conf.option_set( "gen" ) && conf.option_set( "pun" ) ) )
    {
        print_usage(argc, argv);
        std::cout << conf.desc << std::endl;
        return 1;
    }
    if ( conf.setup_configuration() == false )
    {
        std::cout << "The configuration could not been setup." << std::endl;
        return 2;
    }
    
    // Compute all bases.
    MonoComplexZStorageOnly monocomplex( conf.genus, conf.num_punctures, conf.sgn_conv,conf.num_threads, conf.num_remaining_threads );
    
    // Save bases to file
    std::string prefix_basis("./cache/bases/");
    prefix_basis += std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_";
    
    for( auto& it : monocomplex.bases )
    {
        // Store the p-th basis.
        auto& p = it.first;
        save_to_file_bz2<MonoBasis>(it.second, prefix_basis + std::to_string(p));
    }
    
    // Compute all differentials.
    std::string prefix_differentials("./cache/differentials/");
    prefix_differentials += std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_";
    MatrixZDontDiagonalize & differential = monocomplex.matrix_complex.get_current_differential();
    for( auto& it : monocomplex.bases )
    {
        auto& p = it.first;
        size_t num_rows = monocomplex.bases[p].size();
        size_t num_cols = monocomplex.bases[p-1].size();
        initialize_zero_matrix(differential, num_rows, num_cols);
        monocomplex.gen_differential( p );
        save_to_file_bz2<MatrixZDontDiagonalize>( monocomplex.matrix_complex.get_current_differential(), prefix_differentials + std::to_string(p) );
        monocomplex.erase_current_differential();
    }
}
