#ifndef MISC_HPP
#define MISC_HPP

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

/**
 *  @returns true iff the file exists.
 */
bool file_exists( std::string path );

/**
 *  touch given file.
 */
bool touch( std::string path );

/**
 *  @returns true iff the directory exists.
 */
bool directory_exists( std::string path );

/**
 *  @returns true iff the directory could be created.
 */
bool create_directory( std::string path );

/**
 *  create the cache directories.
 */
bool create_working_directories( bool print_status_messages = true );

/**
 * @returns the current date using the POSIX convention. YYYY-MMM-DD HH:MM:SS
 */
std::string current_date();

double current_memory_usage_in_mb();

// Three template functions that are used to print tuples nicely.
template < size_t n, typename... T >
typename std::enable_if< ( n >= sizeof...(T) ) >::type print_tuple( std::ostream&   , const std::tuple< T... >& )
{}

template < size_t n, typename... T >
typename std::enable_if< ( n  < sizeof...(T) ) >::type print_tuple( std::ostream& os, const std::tuple< T... >& tup )
{
    if ( n != 0 )
    {
        os << ", ";
    }
    os << std::get<n>( tup );
    print_tuple< n+1 >( os, tup );
}

template < typename... T >
std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup)
{
    os << "[";
    print_tuple< 0 >( os, tup );
    return os << "]";
}

#endif // MISC_HPP
