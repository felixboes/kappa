#include <future>
#include <fstream>
#include <iostream>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg (-r|-n arg)" << std::endl;
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

    MonoComplexT monocomplex( conf.genus, conf.num_punctures, conf.sgn_conv, conf.num_threads, conf.num_remaining_threads);
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
        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        atomic_uint& current_rank = monocomplex.matrix_complex.get_diagonalizer().current_rank;
        uint32_t max_possible_rank(0);

        // Generate a single differential.
        std::cout << "Constructing the " << p << "-th differential of size " << monocomplex.basis_complex[p].size() << " x " << monocomplex.basis_complex[p-1].size();
        std::cout.flush();
        ofs << "Constructing the " << p << "-th differential of size " << monocomplex.basis_complex[p].size() << " x " << monocomplex.basis_complex[p-1].size();

        measure_duration = Clock();
        monocomplex.gen_differential(p);

        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        
        std::cout << "Matrix: " << std::endl
                  << monocomplex.get_current_differential() << std::endl;

        max_possible_rank = std::min( monocomplex.matrix_complex.num_rows(), monocomplex.matrix_complex.num_cols() );
        if( (uint32_t)homology.get_kern(p-1) > 0 )
        {
            max_possible_rank = std::min( max_possible_rank, (uint32_t)homology.get_kern(p-1) );
        }

        // Compute the induced homology.
        measure_duration = Clock(); // Measure duration.
        atomic_uint state(0);   // Set state to 1 iff kernel and torsion are computed. This is done to terminate the 'monitoring thread'.

        // Diagonalzing thread.
        auto partial_homology_thread = std::async( std::launch::async, [&]() -> HomologyField
        {
            // Always use one thread for diagonalizing at the moment!
            auto ret = monocomplex.matrix_complex.compute_current_kernel_and_torsion( p );
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
        
        
        std::cout << "Matrix: " << std::endl
                  << monocomplex.get_current_differential() << std::endl;
    }

    homology.erase_tors( conf.start_p - 1 );
    homology.erase_kern( conf.end_p + 1 );

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
    std::cout.setf(std::ios::unitbuf);
    // Parse configuration from command line arguments.
    SessionConfig conf(argc, argv);
    if( conf.option_set("help") )
    {
        print_usage(argc, argv);
        std::cout << conf.desc << std::endl;
        return 0;
    }

    if( ! ( conf.option_set( "gen" ) && conf.option_set( "pun" ) && ( conf.option_set( "rat" ) || conf.option_set( "fin" ) ) ) )
    {
        print_usage(argc, argv);
        std::cout << conf.desc << std::endl;
        return 1;
    }

    if ( conf.setup_configuration() == false )
    {
        std::cout << "The configuration could not been setup." << std::endl;
        return 2;
    }

    if (conf.parallel == true)
    {
        Tuple::parallel_case();
    }
    else
    {
        Tuple::radial_case();
    }

    // We may start with the computations.
    if(conf.rational == true)
    {
        compute_homology< MonoComplexQ >( conf, argc, argv );
    }
    else if (conf.prime == 2)
    {
        compute_homology< MonoComplexBool > ( conf, argc, argv);
    }
    else
    {
        compute_homology< MonoComplexZm >( conf, argc, argv );
    }

    return 0;
}
