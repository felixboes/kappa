#ifndef PARALLELIZATION_HPP
#define PARALLELIZATION_HPP

#include <atomic>
#include <future>

typedef std::atomic<uint32_t> atomic_uint;

template< typename Res >
bool thread_running( std::future< Res > & current_thread );

#include "parallelization.ipp"

#endif // PARALLELIZATION_HPP
