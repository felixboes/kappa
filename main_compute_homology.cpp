#include <future>
#include <fstream>
#include <iostream>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures (rational:0 | Z_r: r > 0) optional: first p for which to compute H_p optional:number of threads for paralellization." << std::endl;
}

template< class MonoComplexT >
void compute_homology( SessionConfig conf, int argc, char** argv )
{
    std::ofstream ofs;
    std::string filename = argv[0];
    
    for( int i = 1; i < argc; ++i )
    {
        filename += std::string("_") + std::string(argv[i]);
    }
    
    ofs.open( filename );
    
    std::cout << "-------- Constructing bases --------" << std::endl;
    ofs << "-------- Constructing bases --------" << std::endl;
    
    // Compute all bases.
    Clock measure_duration;
    std::cout << "Constructing bases";
    std::cout.flush();
    ofs << "Constructing bases";
    
    MonoComplexT monocomplex( conf.genus, conf.num_punctures, conf.sgn_conv );
    typename MonoComplexT::HomologyType homology;
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    std::cout << std::endl;
    std::cout.flush();
    ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    ofs << std::endl;
    
    std::cout << "-------- Computing Homology --------" << std::endl;
    ofs << "-------- Computing Homology --------" << std::endl;
    
    // Compute all differentials and homology consecutively.
    for( auto& it : monocomplex.basis_complex )
    {
        int32_t p = it.first;
        if ( p < conf.start_p )
        {
            continue;
        }
        atomic_uint current_rank(0);
        uint32_t max_possible_rank(0);
        
        // Generate a single differential.
        std::cout << "Constructing the " << p << "-th differential";
        std::cout.flush();
        ofs << "Constructing the " << p << "-th differential";
        
        measure_duration = Clock();
        monocomplex.gen_differential(p);
        
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        
        max_possible_rank = std::min( monocomplex.matrix_complex[p].size1(), monocomplex.matrix_complex[p].size2() );
        if( (uint32_t)homology.get_kern(p-1) > 0 )
        {
            max_possible_rank = std::min( max_possible_rank, (uint32_t)homology.get_kern(p-1) );
        }
        
        // Compute the induced homology.
        measure_duration = Clock(); // Measure duration.
        atomic_uint state(0);   // Set state to 1 iff kernel and torsion are computed. This is done to terminate the 'monitoring thread'.

        // Diagonalzing thread.
        auto partial_homology_thread = std::async( std::launch::async, [&]()
        {
            auto ret = monocomplex.matrix_complex.compute_kernel_and_torsion( p, current_rank, conf.num_threads );
            state = 1;
            return ret;
        } );

        // Monitoring thread.
        auto monitor_thread = std::async( std::launch::async, [&]()
        {
            while( state != 1 )
            {
                std::cout << "Diagonalization " << current_rank << "/" << max_possible_rank << "\r";
                std::cout.flush();
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
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
        std::cout << "    dim(H_" << (int32_t)(p-1) << ") = " << (int32_t)(homology.get_kern(p-1) - homology.get_tors(p-1))
                  << "; dim(im D_" << (int32_t)(p) << ") = " << (int32_t)(homology.get_tors(p-1))
                  << "; dim(ker D_" << (int32_t)(p) << ") = " << (int32_t)(homology.get_kern(p)) << std::endl;
        std::cout << std::endl;
        std::cout.flush();
        ofs << "Diagonalization done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << "    dim(H_" << (int32_t)(p-1) << ") = " << (int32_t)(homology.get_kern(p-1) - homology.get_tors(p-1))
            << "; dim(im D_" << (int32_t)(p) << ") = " << (int32_t)(homology.get_tors(p-1))
            << "; dim(ker D_" << (int32_t)(p) << ") = " << (int32_t)(homology.get_kern(p)) << std::endl;
        ofs << std::endl;
        
        // Delete the differential.
        monocomplex.erase_differential(p);
    }
    
    // Print status message.
    std::cout << std::endl;
    std::cout << "------------  Results   ------------" << std::endl;
    std::cout << std::endl;
    ofs << std::endl;
    ofs << "------------  Results   ------------" << std::endl;
    ofs << std::endl;
    
    // Print results.
    std::cout << homology << std::endl;
    ofs << homology << std::endl;
    
    ofs.close();
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
        compute_homology< MonoComplexQ >( conf, argc, argv );
    }
    else
    {
        compute_homology< MonoComplexZm >( conf, argc, argv );
    }
    
    return 0;
}
