#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>
#include <fstream>
#include <string>

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
void create_cache_directories();

#endif // MISC_HPP
