// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#ifndef MISC_HPP
#define MISC_HPP

#include <ctime>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/version.hpp>

#include <bzlib.h>
#include <sys/resource.h>
#include <unistd.h>

#include <libhomology/homology.hpp>

#include "tuple.hpp"

extern const char* program_version_by_git; // will be generated by the make file.

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

/**
 *  Limit the amount of memory that might be allocated.
 *  @warning According to 'man 7 setrlimit',
 *  on 32bit systems the limit is either <= 2GiB or unlimited.
 *  This might result in unwanted behaviour.
 *
 *  int setrlimit(int resource, const struct rlimit *rlim);
 *  
 *  RLIMIT_AS ...
 *  Since the value is a long, on machines with a 32-bit long either this
 *  limit is at most 2 GiB, or this resource is unlimited.
**/
void limit_memory( int32_t percent );

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
 *  @returns the tex code for a given list of cells.
**/
std::string tex_cell( const std::list<Tuple>& cells );

/**
 *  @returns the end of the tex file.
**/
std::string tex_end();
 
/**
 *  @returns program and library versions as well as the date time.
 *  @note In order to knwo the version of the bzip2 library, one needs the header file 'bzlib_private.h' beeing not part of the library package on some distros.
 *        Nontheless, it is part of the source files of the bzip2 package.
 *        The version is given by the macro 'BZ_VERSION'.
**/
std::string kappa_version( int argc = 0, char** argv = nullptr );

template< class CoefficientT >
std::string filename_prefix_differentials( const uint32_t g, const uint32_t m, const bool radial );

template<>
std::string filename_prefix_differentials<Q>( const uint32_t g, const uint32_t m, const bool radial );

template<>
std::string filename_prefix_differentials<Zm>( const uint32_t g, const uint32_t m, const bool radial );

/**
 *  lists all set partitions by declaring to which set an element belongs.
**/
std::vector< std::vector< std::vector< size_t > > > list_set_partitions( size_t n );

/**
 *  lists all connected partitions by declaring to which set an element belongs.
**/
std::vector< std::vector< std::vector< size_t > > > list_connected_partitions( size_t n );

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
