#include <iostream>
#include <string>
#include <vector>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg --first_basis arg --last_basis arg" << std::endl;
}

void print_basis( MonoBasis& M )
{
    std::cout << "Number of basis elements: " << M.size() << std::endl;
    for( auto tup : M.basis )
    {
        Tuple::ConnectedComponents comp = tup.connected_components();
        
        std::cout << " " << tup << "; Listing components: ";
        
        std::vector<bool> visited( comp.size()+1, false );    // visited[0] will not used
        for( uint32_t i = 1; i < comp.size(); )
        {
            std::cout << comp[i] << ":" << i;
            visited[i] = true;
            
            for( int32_t j = i+1; j < comp.size(); ++j )   // mark all visited vertices
            {
                if( comp[j] == comp[i] )
                {
                    visited[j] = true;
                    std::cout << " " << j;
                }
            }
            std::cout << "; ";
    
            // find the next unvisited vertex
            for( ++i; i <= comp.size() && visited[i]; ++i )
            {
            }
        }
        
        std::cout << std::endl;
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
    
    if( ! ( conf.option_set( "gen" ) && conf.option_set( "pun" ) && conf.option_set( "first_basis" ) && conf.option_set( "last_basis" ) ) )
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

    for( auto i = std::max<uint32_t>( conf.first_basis, conf.num_punctures+1 ); i <= std::min<uint32_t>( conf.last_basis, 4*conf.genus + 2*conf.num_punctures ); ++i )
    {
        #ifndef WE_USE_AN_OLD_COMPILER_THAT_DOES_NOT_SUPPORT_ALL_CPP_ELEVEN_FEATURES_OR_OPTIMIZATION
        MonoBasis B = load_from_file_bz2<MonoBasis>("./cache/bases/" + std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_" + std::to_string(i) );
        #else
        MonoBasis B;
        load_from_file_bz2<MonoBasis>( B, "./cache/bases/" + std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_" + std::to_string(i) );
        #endif
    
        print_basis(B);
    }

    return 0;
}
