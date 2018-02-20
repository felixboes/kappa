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


#ifndef CLOCK_HPP
#define CLOCK_HPP

// Description:
//
// This header defines a clock used for time measurements.

#include <chrono>

/**
 *  The class Clock is used for time measurements.
 */ 
class Clock
{
public:
    Clock();    ///< Construcs a Clock and starts the time measurement.
    
    double duration();  ///< @returns the number of secounds since construciton.
protected:
/// @bug Since we want to use c++11 but the g++ version of ubuntu 12.04 does not yet support all features we had to alter the code.
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
    typedef std::chrono::steady_clock std_clock;
#else
    typedef std::chrono::system_clock std_clock;
#endif
    std_clock::time_point measure_duration;

};

#endif // CLOCK_HPP
