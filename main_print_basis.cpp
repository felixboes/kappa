#include <iostream>
#include <string>

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
        std::cout << tup << std::endl;
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
