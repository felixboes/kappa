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


bool create_working_directories( bool print_status_messages )
{
    bool status = true;
    std::list<std::string> dirs = {
        "./cache/",
        "./cache/bases_parallel/",
        "./cache/bases_radial/",
        "./cache/differentials_parallel",
        "./cache/differentials_radial",
        "./cache/list_of_files_that_should_not_be_overwritten",
        "./results"
    };
    
    for( const auto& dir : dirs )
    {
        if( directory_exists( dir ) == false )
        {
            if( print_status_messages == true )
            {
                std::cout << "Creating directory '" << dir << "'.";
            }
            if( create_directory( dir ) == false )
            {
                if( print_status_messages == true )
                {
                    std::cout << " Creation failed." << std::endl;
                }
                status = false;
            }
            else
            {
                std::cout << std::endl;
            }
        }
    }
    
    return status;
}

std::string current_date(){
    return boost::posix_time::to_simple_string( boost::posix_time::ptime(boost::posix_time::second_clock::local_time()) );
}
