#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <boost/serialization/set.hpp>

#include <homology.hpp>

#include "serialization.hpp"
#include "kappa.hpp"
#include "tuple.hpp"

void test1(int argc, char** argv)
{
    if(argc > 2)
    {
        Zm::set_modulus(atoi(argv[1]), 1);
        MatrixZm M( atoi(argv[2]), atoi(argv[2]) );
        M(0,0) = 13;
        
        std::cout << M << std::endl;
        
        save_to_file_bz2( M, "test_out_file" );
        std::cout << load_from_file_bz2<MatrixZm>( "test_out_file" ) << std::endl;
    }
    if(argc > 1)
    {
        int32_t t = atoi(argv[1]);
        
        save_to_file_bz2( t, "test_out_file" );
        load_from_file_bz2<int32_t>( "test_out_file" );
    }
    else
    {
        Tuple t(3);
        t.p=3;
        t.id= 19;
        t[3]=Transposition(3,1);
        t[2]=Transposition(3,2);
        t[1]=Transposition(2,1);
        
        save_to_file_bz2(t, "test_out_file" );
        std::cout << load_from_file_bz2<Tuple>( "test_out_file" ) << std::endl;
    }
}

void test2(int argc, char** argv)
{
    if(argc > 5)
    {
        int32_t g = atoi(argv[1]);
        int32_t m = atoi(argv[2]);
        int32_t prime = atoi(argv[3]);
        int32_t p = atoi(argv[4]);
        
        Zm::set_modulus(prime,1);
        
        if( argv[5][0] == 's' )
        {
            MonoComplexZm M( g, m );
            M.gen_differential(p);
            std::cout << M.matrix_complex[p] << std::endl;
            
            save_to_file_bz2( M.matrix_complex[p], "test_out_file" );
        }
        else
        {
            std::cout << load_from_file_bz2<MatrixZm>( "test_out_file" ) << std::endl;
        }
    }
}

void test3(int argc, char** argv)
{   
    if(argc < 4)
    {
        std::cout << "g m p" << std::endl;
        return;
    }
    std::string t("./cache/bases/");
    t += std::to_string(atoi(argv[1])) + "_" + std::to_string(atoi(argv[2])) + "_" + std::to_string(atoi(argv[3]));
    
    auto ba = load_from_file_bz2<MonoBasis>(t).basis;
    std::cout << ba.size() << std::endl;
    
    for( auto& it : ba )
    {
        std::cout << it << std::endl;
    }
    std::cout << std::endl;    
}

int main(int argc, char** argv)
{
    //test1(argc, argv);
    //test2(argc, argv);
    test3(argc, argv);
    
    return 0;
}
