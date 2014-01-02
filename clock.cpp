#include "clock.hpp"

Clock::Clock()
{
    measure_duration = chrn::steady_clock::now();
}

double Clock::duration()
{
    return chrn::duration<double>(chrn::steady_clock::now() - measure_duration).count();
}
