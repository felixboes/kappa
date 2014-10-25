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

#endif // MISC_HPP