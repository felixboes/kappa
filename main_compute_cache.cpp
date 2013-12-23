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
    SessionConfig conf(argc, argv);
    if( conf.valid == false )
    {
        print_usage(argc, argv);
        return 1;
    }
    
    conf.setup_configuration();
    
    if(conf.rational == true)
    {
        MonoComplexZStorageOnly monocomplex( conf.genus, conf.num_punctures );
        monocomplex.gen_differentials();
        
//        for( auto& it : monocomplex.basis_complex )
//        {
//            auto& p = it.first;
//            std::string prefix_basis("./cache/bases/");
//            std::string prefix_differentials("./cache/differentials/");
            
//            std::string t ("test_out_file_");
//            t += std::to_string(p);
//            save_to_file_bz2<MonoBasis>(it.second, t.c_str());
//            std::cout << it.second.size() << std::endl;
//        }
        
        save_to_file_bz2<MatrixZDontDiagonalize>( monocomplex.matrix_complex[3], "test_out_file" );
        
        return 0;
    }
}
