#include <future>
#include <iostream>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures (rational:0 | Z_r: r > 0)" << std::endl;
}

template< class MonoComplexT >
void compute_homology( SessionConfig conf )
{
    std::cout << "-------- Constructing bases --------" << std::endl;
    
    // Compute all bases.
    MonoComplexT monocomplex( conf.genus, conf.num_punctures );
    typename MonoComplexT::HomologyType homology;
    
    std::cout << std::endl;
    std::cout << "-------- Computing Homology --------" << std::endl;
    
    // Compute all differentials and homology consecutively.
    for( auto& it : monocomplex.basis_complex )
    {
        int32_t p = it.first;
        atomic_uint current_rank(0);
        uint32_t max_possible_rank(0);
        
        // Generate a single differential.
        monocomplex.gen_differential(p);
        max_possible_rank = std::min( monocomplex.matrix_complex[p].size1(), monocomplex.matrix_complex[p].size2() );
        
        // Compute the induced homology.
        Clock measure_duration;
        
        // Diagonalzing thread.
        auto partial_homology_thread = std::async( std::launch::async, [&]()
        {
            return monocomplex.matrix_complex.compute_kernel_and_torsion( p, current_rank );
        } );
        
        // Monitoring thread.
        auto monitor_thread = std::async( std::launch::async, [&]()
        {
            while( thread_running( partial_homology_thread ) )
            {
                std::cout << "Diagonalization " << current_rank << "/" << max_possible_rank << "\r";
                std::cout.flush();
            }
        } );
        
        // Wait for threads to terminate.
        auto partial_homology = partial_homology_thread.get();
        monitor_thread.get();
        
        // Save results.
        homology.set_kern( p, partial_homology.get_kern(p) );
        homology.set_tors( p-1, partial_homology.get_tors(p-1) );
        
        // Print status message.
        std::cout << "Diagonalization done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout << "    dim(H_" << (int32_t)(p-1) << ") = " << (int32_t)(homology.get_kern(p-1) - homology.get_tors(p-1)) << std::endl;
        std::cout << std::endl;
        
        // Delete the differential.
        monocomplex.erase_differential(p);
    }
    
    // Print status message.
    std::cout << std::endl;
    std::cout << "------------  Results   ------------" << std::endl;
    std::cout << std::endl;
    
    // Print results.
    std::cout << homology << std::endl;
}

int main(int argc, char** argv)
{
    // Parse configuration from command line arguments.
    SessionConfig conf(argc, argv);
    if( conf.valid == false )
    {
        print_usage(argc, argv);
        return 1;
    }
    
    if ( conf.setup_configuration() == false )
    {
        std::cout << "The configuration could not been setup." << std::endl;
        return 2;
    }
    
    // We may start with the computations.
    if(conf.rational == true)
    {
        compute_homology< MonoComplexQ >( conf );
    }
    else
    {
        compute_homology< MonoComplexZm >( conf );
    }
    
    return 0;
}
