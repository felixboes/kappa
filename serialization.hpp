#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <iostream>
#include <fstream>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>


/**
 *  Save a given class T to filename.bz2.
 */
template <class T>
void save_to_file_bz2(T& t, std::string filename){
    // Open binary file 'filename.bz2'.
    std::ofstream ofs( (std::string(filename) + ".bz2").c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    
    // Wrap the 'compress with bzip2'-filter around ofs.
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    out.push(boost::iostreams::bzip2_compressor());
    out.push(ofs);
    
    // Open binary_archive in order to serialize the given class.
    boost::archive::binary_oarchive oa(out);
    
    // Put class into the binary_archive.
    oa << t;
}

/**
 *  Save a given class T from filename.bz2.
 */
template <class T>
T load_from_file_bz2(char* filename)
{
    // Open binary file 'filename.bz2'.
    std::ifstream ifs( (std::string(filename) + ".bz2").c_str(), std::ios::in|std::ios::binary);
    
    // Wrap the 'decompress with bzip2'-filter around ifs.
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(ifs);
    
    // Open binary_archive and deserialize the given class.
    T t;
    boost::archive::binary_iarchive(in) >> t;
    
    return t;
}
#endif // SERIALIZATION_HPP
