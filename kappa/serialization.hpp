#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <libhomology/homology.hpp>

/// @warning: unorders_set is not yet supported by boost::serialization.

/**
 *  Save a given class T to filename.bz2.
 */
template <class T>
void save_to_file_bz2( T& t, std::string filename, bool print_duration = true )
{
    if( print_duration == true )
    {
        std::cout << "Saving '" + filename + ".bz2'";
        std::cout.flush();
    }
    Clock measure_duration;
    
    // Open binary file 'filename.bz2'.
    std::ofstream ofs( (filename + ".bz2").c_str(), std::ios::out|std::ios::binary|std::ios::trunc );    
    
    // Wrap the 'compress with bzip2'-filter around ofs.
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    out.push(boost::iostreams::bzip2_compressor());
    out.push(ofs);
    
    // Open binary_archive in order to serialize the given class.
    boost::archive::binary_oarchive oa(out);
    
    // Put class into the binary_archive.
    oa << t;
    
    if( print_duration == true )
    {
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
    }
}

/**
 *  Load a given class T from filename.bz2.
 */
#ifndef WE_USE_AN_OLD_COMPILER_THAT_DOES_NOT_SUPPORT_ALL_CPP_ELEVEN_FEATURES_OR_OPTIMIZATION
template <class T>
T load_from_file_bz2( std::string filename, bool print_duration = true)
#else
template <class T>
void load_from_file_bz2( T& t, std::string filename, bool print_duration = true)
#endif
{
    if( boost::filesystem::exists(filename + ".bz2") && boost::filesystem::is_regular(filename + ".bz2") )
    {
        if( print_duration == true )
        {
            std::cout << "Loading '" + filename + ".bz2'";
            std::cout.flush();
        }
        Clock measure_duration;
        
        // Open binary file 'filename.bz2'.
        std::ifstream ifs( (filename + ".bz2").c_str(), std::ios::in|std::ios::binary );
        
        // Wrap the 'decompress with bzip2'-filter around ifs.
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
        in.push(boost::iostreams::bzip2_decompressor());
        in.push(ifs);
        
        // Open binary_archive and deserialize the given class.
        #ifndef WE_USE_AN_OLD_COMPILER_THAT_DOES_NOT_SUPPORT_ALL_CPP_ELEVEN_FEATURES_OR_OPTIMIZATION
        T t;
        #endif
        boost::archive::binary_iarchive(in) >> t;
        
        if( print_duration == true )
        {
            std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            std::cout.flush();
        }
        #ifndef WE_USE_AN_OLD_COMPILER_THAT_DOES_NOT_SUPPORT_ALL_CPP_ELEVEN_FEATURES_OR_OPTIMIZATION
        return t;
        #endif
    }
    else
    {
        std::cout << "This '" << filename + ".bz2' is not a file." << std::endl;
        std::cout.flush();
        #ifndef WE_USE_AN_OLD_COMPILER_THAT_DOES_NOT_SUPPORT_ALL_CPP_ELEVEN_FEATURES_OR_OPTIMIZATION
        return T();
        #endif
    }
}
#endif // SERIALIZATION_HPP
