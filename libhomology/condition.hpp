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


#ifndef CONDITION_HPP
#define CONDITION_HPP

#include <condition_variable>
#include <mutex>

/*! \class Condition
 *  \brief Basic boolean condition that can be signalled and waiter for
 *
 *  A \c Condition models a boolean condition wrapped around a single
 *  \c std::condition_variable. A condition can be signalled by calling
 *  \c signal telling waiting idle threads the condition is satisfied. A thread
 *  can call \c wait in order to stay idle until the condition is met (or more
 *  precisely until the thread is notified that the condition is met).
 *
 *  If signalling happens before any threads starts waiting, the above mechanism
 *  still works: A \c Condition stores that is has been signalled until a thread
 *  receives the notification by calling \c wait. Signalling more than once
 *  before the signal is waited for will have the same effect as signalling
 *  once.
 *
 *  All changes on the modelled condition are performed locked in order to
 *  facilitate a safe concurrent use of \c wait and \c signal.
 *
 *  Here is an exemplary usage adapted from
 *  <a href="http://en.cppreference.com/w/cpp/thread/condition_variable">
 *  http://www.cppreference.com</a>:
 *  \code
 *  std::string  data;
 *  Condition    wait_cond;
 *  Condition    done_cond;
 *
 *  void worker_func()
 *  {
 *      // Wait until main() sends data
 *      wait_cond.wait();
 *      data += " after processing";
 *
 *      // Send data back to main()
 *      std::cout << "Worker thread signals data processing completed\n";
 *      done_cond.signal();
 *
 *      // Unlocking was done by the Condition variables. Nothing further here!
 *  }
 *
 *  int main()
 *  {
 *      // Start the function worker_thread in a parallel thread
 *      std::thread worker(worker_thread);
 *
 *      // Signal the worker_func to actually start working
 *      data = "Example data";
 *      std::cout << "main() signals data ready for processing\n";
 *      wait_cond.signal();
 *
 *      // Wait for the worker to finish
 *      done_cond.wait();
 *      std::cout << "Back in main(), data = " << data << '\n';
 *
 *      worker.join();
 *  }
 *  \endcode
**/
class Condition
{
public:
    //! Basic constructor initiallizing non-signalled
    Condition();

    /*! \brief Wait for the condition to be signalled.
     *
     *  Wait until the condition is met (or more precisely until the calling
     *  thread is signalled that the condition is met). If the \c Condition is
     *  not signalled anymore, any waiting threads will remain waiting
     *  infinitely. After receiving the signal, the \c Condition stores that
     *  signalling is received (ie. sets \c triggered to \c false).
    **/
    void wait();

    /*! \brief Notify one waiting thread about the condition being met
     *
     *  Notify one waiting thread that the condition is met now. If multiple
     *  threads are waiting, only one of them will wake up.
     *
     *  The fact that signalling took place is stored (ie. \c triggered is set
     *  to \c true). This allows other threads to receive the signal even if
     *  their \c wait call succeeds this \c signal call (ie. a thread intended
     *  to wait for \c signal invokes \c wait after \c signal is performed).
    **/
    void signal();

protected:
    /*! Mutex for locking whenever condition changes **/
    std::mutex                    mtx;
    //! Condition Variable for waiting 'till condition is signalled
    std::condition_variable       cv;
    //! Boolean to distinguish between user signalling and spurious wake
    bool                          triggered;
}; // class Condition


#endif // CONDITION_HPP
