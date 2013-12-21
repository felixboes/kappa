/*  Minimum working example.
 *  libbzip2-dev must be installed.
 *  Complile with
 *      g++ -std=c++11 test_serialization.cpp -lboost_iostreams -lboost_serialization
 *  Run with
 *      ./a.out number_to_save
 */
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <iostream>
#include <fstream>
#include <string>

void save(int32_t o)
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
    oa << o;
}

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
    
    int i;
    ia >> i;
    
    // Output result.
    std::cout << i << std::endl;
}

int main(int argc, char** argv)
{
    if(argc > 1)
    {
        save( atoi(argv[1]) );
    }
    
    load();
    return 0;
}
