#include <iostream>
#include <string>
#include <vector>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures num_module" << std::endl;
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
    if(argc < 4)
    {
        print_usage(argc, argv);
        return 1;
    }

    MonoBasis B = load_from_file_bz2<MonoBasis>("./cache/bases/" + std::string(argv[1]) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) );

    print_basis(B);

    return 0;
}
