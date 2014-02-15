template< typename Res >
bool thread_running( std::future< Res > & thread_to_check )
{
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
    // Check wether the thread is valid. If so wait at most 100 milliseconds and check why the functions returned. Here std::future_status::timeout means that it is still working.
    return thread_to_check.valid() && thread_to_check.wait_for( std::chrono::milliseconds(100) ) == std::future_status::timeout;
#elif __GNUC__ == 4 && __GNUC_MINOR__ > 5
    // http://gcc.gnu.org/onlinedocs/gcc-4.6.4/libstdc++/api/a00488.html
    return thread_to_check.valid() &&  ! thread_to_check.wait_for( std::chrono::milliseconds(100) );
#endif
}
