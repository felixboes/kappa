#ifndef PARALLELIZATION_HPP
#define PARALLELIZATION_HPP

// Description:
//
// This header defines features can that are used in multithreaded functions.

#include <atomic>
#include <future>

/**
 *  In order to work with thread save integers we define an atomic uint.
 */
typedef std::atomic<uint32_t> atomic_uint;

/**
 *  In order to work with several g++-versions that do not yet support full c++11 namely g++ <= 4.6 on ubuntu 12.04 or c++ > 4.6,
 *  you can use this function to see wether a thread is still running or not.
 *  @returns true iff the thread is still running.
 */
template< typename Res >
bool thread_running( std::future< Res > & current_thread );

#include "parallelization.ipp"

#endif // PARALLELIZATION_HPP
