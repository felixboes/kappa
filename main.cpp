#include <iostream>

#include <homology.hpp>

#include "kappa.hpp"

struct session_config
{
    session_config() : genus(0), num_punctures(0), rational(0), prime(2), valid(false) {}
    
    operator bool() {return valid;}
    
    uint32_t genus;
    uint32_t num_punctures;
    bool rational;
    uint32_t prime;
    bool valid;
};

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures (rational:0 | Z_r: r > 0)" << std::endl;
}

session_config process_command_line_arguments(int argc, char** argv)
{
    if( argc < 4)
    {
        return session_config();
    }
    session_config conf;
    
    conf.genus = atoi(argv[1]);
    conf.num_punctures = atoi(argv[2]);
    conf.rational = (atoi(argv[3]) == 0 ? true : false);
    conf.prime = atoi(argv[3]);
    conf.valid = true;
    
    return conf;
}

bool setup_configuration(session_config conf)
{
    if( conf )
    {
        if( conf.rational == false )
        {
            Zm::set_modulus(conf.prime, 1);
        }
        return true;
    }
    else
    {
        return false;
    }
}

int main(int argc, char** argv)
{
    session_config conf = process_command_line_arguments(argc, argv);
    if( setup_configuration(conf) == false )
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
