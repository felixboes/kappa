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


template< typename Res >
bool thread_running( std::future< Res > & thread_to_check )
{
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
    // Check wether the thread is valid. If so wait at most 100 milliseconds and check why the functions returned. Here std::future_status::timeout means that it is still working.
    return thread_to_check.valid() && thread_to_check.wait_for( std::chrono::milliseconds(100) ) == std::future_status::timeout;
#elif __GNUC__ == 4 && __GNUC_MINOR__ > 5
    // http://gcc.gnu.org/onlinedocs/gcc-4.6.4/libstdc++/api/a00488.html
    return thread_to_check.valid() &&  ! thread_to_check.wait_for( std::chrono::milliseconds(100) );
#endif
}
