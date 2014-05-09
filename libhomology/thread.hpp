#ifndef THREAD_HPP
#define THREAD_HPP

#include <functional>


/*! \class Thread
 *  \brief Basic Thread class allowing Lambda executions in concurrency
 *
 *  This class differs only sightly from \c std::thread. The latter alike, it
 *  can be started immediately upon creation (using an appropriate \c StartFlag)
 *  and waited for to finish working. But in contrast to \c std::thread it can
 *  be restarted (on different input or with a different task).
 *
 *  This \c Thread maintains an \c std::thread internally that will turn idle
 *  when nothing has to be performed. If told so, the internal \c std::thread
 *  will start executing a given \c Lambda in parallel, turning idle afterwards
 *  again.
 *
 *  \warning It is the users responsibility to ensure that each execution of
 *  work is finished, before changing the task or re-executing. Neglecting this
 *  leads to undefined behaviour.
 *
 *  Let's see an example:
 *  \code
 *  Thread::Lambda work = []{std::this_thread::sleep_for(std::chrono::seconds(10));};
 *  Thread worker(work); // Creates an idle Thread
 *
 *  worker.execute(); //  Start executing the given Lambda work in parallel
 *
 *  // Next two lines are forbidden; User's have to ensure this situation does
 *  // not occur before the work is completely done
 *  worker.execute();
 *  worker.set_function([]{std::this_thread::sleep_for(std::chrono::seconds(10));});
 *
 *  worker.wait(); // Terminates when the execute call finished
 *  \endcode
 *
 *  So concluding as a rule of thumb, any call of \c execute should be followed
 *  by \c wait before any other member function is invoked (while of course
 *  other threads may do other work).
 *
 *  \note Implementation uses the PIMPL design idiom (in c++11 style)
**/
class Thread
{
public:
    /*! \brief Basic Lambda typedef for all executions to be performed
     *
     *  Every work that has to be done must be passed as a \c Lambda. This
     *  simplifies the interface.
    **/
    typedef std::function<void ()> Lambda;

    //! Basic enum to increase readability when calling the constructor
    enum StartFlag
    {
        start_immediately, //!< Start executing the lambda immediately
        postpone_start     //!< Postpone lambda invocation 'till execute is called
    }; // enum StartFlag

    /*! \brief Basic Constructor initializing members.
     *
     *  If <tt>start == start_immediately</tt>, \c f will be invoked
     *  immediately.
    **/
    Thread(Lambda const & f, StartFlag start = postpone_start);

    //! Move constructor in order to maintain Threads in STL containers.
    Thread(Thread && thread);

    /*! \brief Waits for the last execution to finish and joins the worker
     *  thread.
     *
     * Note: While being implemented trivially, this method must not be compiler
     * generated (as otherwise unique_ptr can't call the destructor as
     * \c Thread::Impl is not yet defined). User defining it after Thread::Impl
     * has been defined works though.
    **/
    ~Thread();

    //! Sets the stored Function
    void set_function(Lambda const &);

    //! Executes stored lamba in parallel
    void execute();

    //! Convenience call \c set_function() immediately followed by \c execute().
    void execute(Lambda const &);

    //! Returns when the function is exited. The thread will turn idle afterwards.
    void wait();

private:
    class Impl;
    std::unique_ptr<Impl>      pimpl;
}; // class Thread


#endif // THREAD_HPP