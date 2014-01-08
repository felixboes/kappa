#ifndef CLOCK_HPP
#define CLOCK_HPP

#include <chrono>

class Clock
{
public:
    Clock();
    
    double duration();
private:
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
    typedef std::chrono::steady_clock std_clock;
#else
    typedef std::chrono::system_clock std_clock;
#endif
    std_clock::time_point measure_duration;

};

#endif // CLOCK_HPP
