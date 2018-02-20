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
 *  In order to work with several g++-versions that do not yet support full c++11 namely g++ <= 4.6 on ubuntu 12.04 or g++ > 4.6,
 *  you can use this function to see wether a thread is still running or not.
 *  @returns true iff the thread is still running.
 *  @warning: The standard claims, that wait_for does return in a defined behaviour if the future type is not valid.
 *  Therefore one would have to guarantee that checking validness and then calling wait_for is atomic.
 *  Better use a workaround with states stored in an atomic_uint.
 */
//template< typename Res >
//bool thread_running( std::future< Res > & current_thread );

//#include "parallelization.ipp"

#endif // PARALLELIZATION_HPP
