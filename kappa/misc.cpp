#include "misc.hpp"

bool file_exists( std::string path )
{
    boost::filesystem::path p(path);
    try
    {
        if( boost::filesystem::exists(p) && boost::filesystem::is_regular_file(p) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch( const boost::filesystem::filesystem_error& ex )
    {
        std::cout << ex.what() << std::endl;
    }
    return false;
}

bool touch( std::string path )
{
    // convert unix path to native path.
    boost::filesystem::path p(path);
    std::string native_path = p.native();
    
    std::ofstream ofs;
    ofs.open(native_path);
    if( ofs.is_open() == false )
    {
        ofs.close();
        return false;
    }
    ofs.close();
    return true;
}

bool directory_exists( std::string path )
{
    boost::filesystem::path p(path);
    try
    {
        if( boost::filesystem::exists(p) && boost::filesystem::is_directory(p) )
        {
           return true;
        }
        else
        {
            return false;
        }
    }    
    catch( const boost::filesystem::filesystem_error& ex )
    {
        std::cout << ex.what() << std::endl;
    }
    return false;
}

bool create_directory( std::string path )
{
    boost::filesystem::path p(path);
    try
    {
        if( boost::filesystem::create_directory(p) == true )
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch( const boost::filesystem::filesystem_error& ex )
    {
        std::cout << ex.what() << std::endl;
    }
    return false;
}


void create_cache_directories()
{
    std::list<std::string> dirs = {
        "./cache/",
        "./cache/bases/",
        "./cache/differentials_parallel",
        "./cache/differentials_radial",
        "./cache/list_of_files_that_should_not_be_overwritten"
    };
    for( const auto& dir : dirs )
    {
        if( directory_exists( dir ) == false )
        {
            create_directory( dir );
        }
    }
}
