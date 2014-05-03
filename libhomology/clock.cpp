#include "clock.hpp"

Clock::Clock()
{
    measure_duration = std_clock::now();
}

double Clock::duration()
{
    return std::chrono::duration<double>(std_clock::now() - measure_duration).count();
}
