#include <future>
#include <fstream>
#include <iostream>

#include <homology.hpp>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg (-r|-n arg)" << std::endl;
}

template< class ClusterSpectralSequenceT >
void compute_css( SessionConfig conf, int argc, char** argv )
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
    
    ClusterSpectralSequenceT cluster_spectral_sequence( conf.genus, conf.num_punctures, conf.sgn_conv );
    typename ClusterSpectralSequenceT::CSSHomologyType homology;
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    std::cout << std::endl;
    std::cout.flush();
    ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    ofs << std::endl;

    std::cout << "-------- Computing E^1-term --------" << std::endl;
    ofs << "-------- Computing E^1-term --------" << std::endl;
    
    // Compute all differentials and homology consecutively.
    for( auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        auto p = basis_it.first;
        auto l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        atomic_uint current_rank(0);
        uint32_t max_possible_rank(0);
        
        // Generate a single differential.
        std::cout << "Constructing the differentials of E^0_{" << p << ",*}.";
        std::cout.flush();
        ofs << "Constructing the differentials of E^0_{" << p << ",*}.";

        measure_duration = Clock();
        cluster_spectral_sequence.gen_differentials(p);
        
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;

        for( auto& l_basis_it : l_bases )
        {
            auto l = l_basis_it.first;

            max_possible_rank = std::min( cluster_spectral_sequence.css_page[p][l].size1(), cluster_spectral_sequence.css_page[p][l].size2() );
            if( (uint32_t)homology[p-1].get_kern(l) > 0 )
            {
                max_possible_rank = std::min( max_possible_rank, (uint32_t)homology[p-1].get_kern(l) );
            }

            // Compute the induced homology.
            measure_duration = Clock(); // Measure duration.
            atomic_uint state(0);   // Set state to 1 iff kernel and torsion are computed. This is done to terminate the 'monitoring thread'.

            // Diagonalzing thread.
            auto partial_homology_thread = std::async( std::launch::async, [&]()
            {
                auto ret = cluster_spectral_sequence.css_page[p].compute_kernel_and_torsion( l, current_rank, conf.num_threads );
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
            homology[p].set_kern( l, partial_homology.get_kern(l) );
            homology[p-1].set_tors( l, partial_homology.get_tors(l-1) ); // Observe that the computed torsion is located in position l-1.

            // Print status message.
            std::cout << "Diagonalization done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            std::cout << "    dim(E^1_{" << (int32_t)(p-1) << "," << (int32_t)(l) << "}) = " << (int32_t)(homology[p-1].get_kern(l) - homology[p-1].get_tors(l))
                      << "; dim(im d^1_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << (int32_t)(homology[p-1].get_tors(l))
                      << "; dim(ker d^1_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << (int32_t)(homology[p].get_kern(l)) << std::endl;
            std::cout.flush();
            ofs << "Diagonalization done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            ofs << "    dim(E_" << (int32_t)(p-1) << "," << (int32_t)(l) << ") = " << (int32_t)(homology[p-1].get_kern(l) - homology[p-1].get_tors(l))
                      << "; dim(im d_" << (int32_t)(p) << "," << (int32_t)(l) << ") = " << (int32_t)(homology[p-1].get_tors(l))
                      << "; dim(ker d_" << (int32_t)(p) << "," << (int32_t)(l) << ") = " << (int32_t)(homology[p].get_kern(l)) << std::endl;
        }
        
        // Delete the differential.
        cluster_spectral_sequence.erase_differentials(p);
    }

    
    // Print status message.
    std::cout << std::endl;
    std::cout << "------------  Results   ------------" << std::endl;
    std::cout << std::endl;
    ofs << std::endl;
    ofs << "------------  Results   ------------" << std::endl;
    ofs << std::endl;
    
    // Print Dimensions of the E^0-page
    std::cout << "Dimensions of the E^0 page" << std::endl;
    std::cout.flush();
    ofs << "Dimensions of the E^0 page" << std::endl;
    for( auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        auto p = basis_it.first;
        auto l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        for( auto& l_basis_it : l_bases )
        {
            auto l = l_basis_it.first;
            auto cur_basis = l_basis_it.second;
            std::cout << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") = " << std::setw(4) << cur_basis.size() << ";    ";
            std::cout.flush();
            ofs << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") = " << std::setw(4) << cur_basis.size() << ";    ";
        }
        std::cout << std::endl;
        std::cout.flush();
        ofs << std::endl;
    }
    
    // Print Dimensions of the E^1-page
    std::cout << std::endl;
    std::cout << "Dimensions of the E^1 page" << std::endl;
    std::cout.flush();
    ofs << std::endl;
    ofs << "Dimensions of the E^1 page" << std::endl;
    for( auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        auto p = basis_it.first;
        auto l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        for( auto& l_basis_it : l_bases )
        {
            auto l = l_basis_it.first;
            std::cout << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") = " << std::setw(4) << (int32_t)(homology[p].get_kern(l) - homology[p].get_tors(l)) << ";    ";
            std::cout.flush();
            ofs << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") = " << std::setw(4) << (int32_t)(homology[p].get_kern(l) - homology[p].get_tors(l)) << ";    ";
        }
        std::cout << std::endl;
        std::cout.flush();
        ofs << std::endl;
    }
    
    
    ofs.close();
}

int main(int argc, char** argv)
{
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
    
    // We may start with the computations.
    if(conf.rational == true)
    {
        compute_css< ClusterSpectralSequenceQ >( conf, argc, argv );
    }
    else
    {
        compute_css< ClusterSpectralSequenceZm >( conf, argc, argv );
    }
    
    return 0;
}

