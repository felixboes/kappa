/*  Minimum working example.
 *  libbzip2-dev must be installed.
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

#include "tuple.hpp"

void save_int(int32_t t)
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

void save_tuple(Tuple t)
{
    // Open binary file 'test_out_file.bz2'.
    std::ofstream ofs("test_out_file.bz2", std::ios::out|std::ios::binary|std::ios::trunc);
    
    // Wrap the 'compress with bzip2'-filter around ofs.
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    out.push(boost::iostreams::bzip2_compressor());
    out.push(ofs);
    
    // Open binary_archive in order to serialize the given tuple.
    boost::archive::binary_oarchive oa(out);
    
    // Put tuple into the binary_archive.
    oa << t;
}

void load_int()
{
    // Open binary file 'test_out_file.bz2'.
    std::ifstream ifs("test_out_file.bz2", std::ios::in|std::ios::binary);
    
    // Wrap the 'decompress with bzip2'-filter around ifs.
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(ifs);
    
    // Open binary_archive in order to deserialize the given integer.
    boost::archive::binary_iarchive ia(in);
    
    int i;
    ia >> i;
    
    // Output result.
    std::cout << i << std::endl;
}

void load_tuple()
{
    // Open binary file 'test_out_file.bz2'.
    std::ifstream ifs("test_out_file.bz2", std::ios::in|std::ios::binary);
    
    // Wrap the 'decompress with bzip2'-filter around ifs.
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(ifs);
    
    // Open binary_archive in order to deserialize the given tuple.
    boost::archive::binary_iarchive ia(in);
    
    Tuple t;
    ia >> t;
    
    // Output result.
    std::cout << t << std::endl;
}

int main(int argc, char** argv)
{
    if(argc > 1)
    {
        save_int( atoi(argv[1]) );
        load_int();
    }
    else
    {
        Tuple t(3);
        t.p=3;
        t.id= 19;
        t[3]=Transposition(3,1);
        t[2]=Transposition(3,2);
        t[1]=Transposition(2,1);
        
        save_tuple(t);
        load_tuple();
    }
    return 0;
}
