#include "condition.hpp"


Condition::Condition()
:
    triggered(false)
{
    // intentionally do nothing!
}

void Condition::wait()
{
    // Note: While waiting, lk is unlocked s.t. signalling can take place
    std::unique_lock<std::mutex> lk(mtx);
    cv.wait(lk, [&]{return triggered;} );
    triggered = false;
    lk.unlock();
}

void Condition::signal()
{
    std::lock_guard<std::mutex> lk(mtx);
    triggered = true;
    cv.notify_one();
}
