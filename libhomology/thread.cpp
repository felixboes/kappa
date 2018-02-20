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


#include <thread>

#include "condition.hpp"
#include "thread.hpp"

static void thread_function(
    Condition &               start,
    Condition &               done,
    bool &                    finish,
    Thread::Lambda const &    func)
{
    while (not finish)
    {
        // Wait until the function should be called. Spurious wakups are handled
        // internally.
        start.wait();

        if (finish)
        {
            return;
        }
        // Actually invoke the function
        func();

        // Show that the function is terminated
        done.signal();
    }
}


////////////////////////////////////////////////////////////////////////////////


class Thread::Impl
{
public:
    Impl(
        Thread::Lambda const &   f,
        Thread::StartFlag        start = Thread::postpone_start);

    ~Impl();

    void execute();

    void wait();

    void set_function(Thread::Lambda const &);

private:
    Thread::Lambda           func;

    Condition                start_execution;
    Condition                execution_done;
    bool                     die;

    std::thread              worker;
}; // class Impl


Thread::Impl::Impl(Thread::Lambda const & f, Thread::StartFlag start)
:
    func(f),
    die(false),
    worker([&]{ thread_function(start_execution, execution_done, die, func); })
{
    if (start == Thread::start_immediately)
    {
        execute();
    }
}

Thread::Impl::~Impl()
{
    die = true;
    start_execution.signal();
    worker.join();
}

void Thread::Impl::execute()
{
    start_execution.signal();
}

void Thread::Impl::wait()
{
    execution_done.wait();
}

void Thread::Impl::set_function(Thread::Lambda const & f)
{
    func = f;
}


//////////////////////////////////////////////////////////////////////////////


Thread::Thread(Lambda const & f, StartFlag start)
:
    pimpl(new Impl(f, start))
{
    // intentionally do nothing else!
}

Thread::Thread(Thread && other)
{
    pimpl = std::move(other.pimpl);
    other.pimpl = nullptr;
}

Thread::~Thread()
{
    // intentionally do nothing!
}

void Thread::set_function(Lambda const & f)
{
    pimpl->set_function(f);
}

void Thread::execute()
{
    pimpl->execute();
}

void Thread::execute(Lambda const & f)
{
    set_function(f);
    execute();
}

void Thread::wait()
{
    pimpl->wait();
}
