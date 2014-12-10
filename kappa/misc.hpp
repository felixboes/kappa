#ifndef MISC_HPP
#define MISC_HPP

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

#include "tuple.hpp"

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

/**
 *  @returns the current memory usage in mb.
 */ 
double current_memory_usage_in_mb();


// Functions to print cells.
/**
 *  @returns the preamble of the tex file.
**/
std::string tex_preamble();

/**
 *  @returns the tex code for a given cell.
**/
std::string tex_cell( const Tuple& cell );

/**
 *  @returns the end of the tex file.
**/
std::string tex_end();
 

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

// Used to print diagonals.
template< typename T >
std::ostream& operator<<( std::ostream& os, const std::list< std::pair<T,T> >& list )
{
    os << "[ ";
    for( const auto& it : list )
    {
        os << "[" << it.first << "," << it.second << "]; ";
    }
    return os << "]";
}

#endif // MISC_HPP
