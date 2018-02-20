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


#include "condition.hpp"

Condition::Condition() :
    triggered(false)
{
    // intentionally do nothing!
}

void Condition::wait()
{
    // Note: While waiting, lk is unlocked s.t. signalling can take place
    std::unique_lock< std::mutex > lk(mtx);   // A unique lock is used since it is needed for condition_variable.
    cv.wait(lk, [&]{return triggered;} );
    triggered = false;
    lk.unlock();
}

void Condition::signal()
{
    std::lock_guard< std::mutex > lk(mtx);
    triggered = true;
    cv.notify_one();
}
