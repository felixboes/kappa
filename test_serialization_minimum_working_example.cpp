/*  Minimum working example.
 *  libbz2-dev must be installed.
 *  Complile with
 *      g++ -std=c++11 test_serialization.cpp tuple.cpp -lboost_iostreams -lboost_serialization
 *  Run with
 *      ./a.out
 *      ./a.out save_a_single_integer
 */
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <iostream>
#include <fstream>

#include <homology.hpp>

#include "kappa.hpp"
#include "tuple.hpp"

template <class T>
void save(T& t)
{
    // Open binary file 'test_out_file.bz2'.
    std::ofstream ofs("test_out_file.bz2", std::ios::out|std::ios::binary|std::ios::trunc);
    
    // Wrap the 'compress with bzip2'-filter around ofs.
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    out.push(boost::iostreams::bzip2_compressor());
    out.push(ofs);
    
    // Open binary_archive in order to serialize the given integer.
    boost::archive::binary_oarchive oa(out);
    
    // Put integer into the binary_archive.
    oa << t;
}

template <class T>
void load()
{
    // Open binary file 'test_out_file.bz2'.
    std::ifstream ifs("test_out_file.bz2", std::ios::in|std::ios::binary);
    
    // Wrap the 'decompress with bzip2'-filter around ifs.
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(ifs);
    
    // Open binary_archive in order to deserialize the given integer.
    boost::archive::binary_iarchive ia(in);
    
    T i;
    ia >> i;
    
    // Output result.
    std::cout << i << std::endl;
}

void test1(int argc, char** argv)
{
    if(argc > 2)
    {
        Zm::set_modulus(atoi(argv[1]), 1);
        MatrixZm M( atoi(argv[2]), atoi(argv[2]) );
        M(0,0) = 13;
        
        std::cout << M << std::endl;
        
        save<MatrixZm> ( M );
        load<MatrixZm> ();
    }
    if(argc > 1)
    {
        int32_t t = atoi(argv[1]);
        save<int32_t>( t );
        load<int32_t>();
    }
    else
    {
        Tuple t(3);
        t.p=3;
        t.id= 19;
        t[3]=Transposition(3,1);
        t[2]=Transposition(3,2);
        t[1]=Transposition(2,1);
        
        save<Tuple>(t);
        load<Tuple>();
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
            
            save<MatrixZm>( M.matrix_complex[p] );
        }
        else
        {
            load<MatrixZm>();
        }
    }
}

int main(int argc, char** argv)
{
    //test1(argc, argv);
    test2(argc, argv);
    return 0;
}
