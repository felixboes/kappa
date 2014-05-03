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
private:
// Since we want to use c++11 but the g++ version of ubuntu 12.04 does not yet support all features we had to alter the code.
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
    typedef std::chrono::steady_clock std_clock;
#else
    typedef std::chrono::system_clock std_clock;
#endif
    std_clock::time_point measure_duration;

};

#endif // CLOCK_HPP
