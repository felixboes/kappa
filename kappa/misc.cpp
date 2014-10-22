#include "misc.hpp"

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
    if( directory_exists("./cache/") == false )
    {
        create_directory("./cache/");
    }
    if( directory_exists("./cache/differentials_parallel") == false )
    {
        create_directory("./cache/differentials_parallel");
    }
    if( directory_exists("./cache/differentials_radial") == false )
    {
        create_directory("./cache/differentials_radial");
    }
    if( directory_exists("./cache/bases") == false )
    {
        create_directory("./cache/bases");
    }
}
