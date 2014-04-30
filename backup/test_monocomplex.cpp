#include <atomic>
#include <chrono>
#include <future>
#include <iomanip>
#include <iostream>

#include <homology.hpp>
#include "monocomplex.hpp"

template< class MatrixType >
bool equals( const MatrixType& A, const MatrixType& B )
{
    if( A.size1() == B.size1() && A.size2() == B.size2() )
    {
        uint32_t rows = A.size1();
        uint32_t cols = A.size2();
        
        for( uint32_t i = 0; i < rows; ++i )
        {
            for( uint32_t j = 0; j < cols; ++j )
            {
                if( A(i,j) != B(i,j) )
                {
                    return false;
                }
            }
        }
        
        return true;
    } 
    
    return false;
}

void test1( int argc, char** argv )
{
    if( argc < 4 )
    {
        std::cout << "Usage: " << argv[0] << " g m differential_to_print (show_basis = 0 or 1) (check_dd_zero = 0 or 1)" << std::endl;
        return;
    }
    
    typedef MonoComplex<ChainComplexQ> MonoComplexQ;
    MonoComplexQ mc( atoi(argv[1]), atoi(argv[2]) ); // p = 2( 1*2+3 )= 10
    
    if( argc >= 5 && atoi(argv[4]) != 0 )
    {
        mc.show_basis( atoi(argv[3]) - 1 );
        mc.show_basis( atoi(argv[3]) );
    }
    
    mc.gen_differential( atoi(argv[3]) );
    
    std::cout << mc.matrix_complex[atoi(argv[3])] << std::endl;
    
    if( argc >= 6 && atoi(argv[5]) != 0 )
    {
        mc.gen_differential( atoi(argv[3]) - 1 );
        std::cout << mc.matrix_complex[atoi(argv[3]) - 1] << std::endl;
        
        std::cout << "Check if dd = 0:" << std::endl;
        std::cout << prod( mc.matrix_complex[atoi(argv[3]) - 1], mc.matrix_complex[atoi(argv[3])] ) << std::endl; 
    }
}

void test2( int argc, char** argv )
{
    if( argc < 4 )
    {
        std::cout << "Usage: " << argv[0] << " g m differential_to_print (show_basis = 0 or 1) (check_dd_zero = 0 or 1)" << std::endl;
        return;
    }
    
    Zm::set_modulus(2,1);
    
    typedef MonoComplex<ChainComplexZm> MonoComplexZm;
    MonoComplexZm mc( atoi(argv[1]), atoi(argv[2]) ); // p = 2( 1*2+3 )= 10
    
    mc.gen_differential( atoi(argv[3]) );
    mc.gen_differential( atoi(argv[3]) + 1 );
    
    atomic_uint current_rank_0(0);
    atomic_uint current_rank_1(0);
    
    // Start a thread that diagonalizes the matrix given by the 3rd command line argument.
    // The process (rows that are diagonalized) is stores in current_rank_0
    // In detail:
    // We call std::async with the 'policy' std::launch::async i.e. we start a new thread.
    // This thread calls our lambda function.
    // The lambda function accesses the variables by reference as it is defined as [&].
    // The return value of std::async is a std::future<void>
    auto diag_thread_0 = std::async(std::launch::async, [&]()
        {
            std::cout << "Diagonalizer 0 initialized." << std::endl;
            std::cout.flush();
            DiagonalizerZm diago;
            diago( mc.matrix_complex[ atoi(argv[3]) ], current_rank_0 );
        } );
    // As above.
    auto diag_thread_1 = std::async(std::launch::async, [&]()
        {
           std::cout << "Diagonalizer 1 initialized." << std::endl;
           std::cout.flush();
           DiagonalizerZm diago;
           diago( mc.matrix_complex[ atoi(argv[3]) + 1 ], current_rank_1 );
        } );
    // Start a thread that prints out the current status.
    auto monitor_thread = std::async(std::launch::async, [&]()
        {
            std::cout << "Monitor initialized." << std::endl;
            std::cout.flush();
            
            // We print the status as long as one of the diagoanlizer threads is working.
            // Therefore we have to check if the std::future(s) are still valid and
            // wait max. 100 second to check wether the thread is still working.
            // Observe: In the printing command, we HAVE to check wether diag_thread_x is valid since
            // after returning from the get below, the behaviour of wait_for() is undefined.
            // See also [Stroustrup page 1242].
            while( thread_running(diag_thread_0) )
            {
                std::cout << "Thread 0:" << std::setw(5) << current_rank_0 << " " <<
                             "Thread 1:" << std::setw(5) << ( thread_running(diag_thread_1) == false ? "done" : std::to_string(current_rank_1.load()) ) << "\r";
                std::cout.flush();
            }
            while( thread_running( diag_thread_1 ) )
            {
                std::cout << "Thread 0:" << std::setw(5) << ( thread_running( diag_thread_0 ) == false ? "done" : std::to_string(current_rank_0.load()) ) << " " <<
                             "Thread 1:" << std::setw(5) << current_rank_1 << "\r";
                std::cout.flush();
            }
        } );

    
    // Wait for thread terminations.
    diag_thread_0.get();
    diag_thread_1.get();
    monitor_thread.get();
    
    std::cout << "Thread 0:" << std::setw(5) << "done" << " " <<
                 "Thread 1:" << std::setw(5) << "done" << std::endl;
}

int main(int argc, char** argv)
{
//    test1(argc, argv); // Test dd = 0
    test2(argc, argv); // Test paralellization.
    return 0;
}
