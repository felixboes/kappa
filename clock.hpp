#ifndef CLOCK_HPP
#define CLOCK_HPP

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
#include <chrono>
namespace chrn = std::chrono;
#else
#include <boost/chrono.hpp>
namespace chrn = boost::chrono;
#endif

class Clock
{
public:
    Clock();
    
    double duration();
private:
    chrn::steady_clock::time_point measure_duration;
};

#endif // CLOCK_HPP
