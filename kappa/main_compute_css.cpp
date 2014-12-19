#include <future>
#include <fstream>
#include <iostream>

#include <libhomology/homology.hpp>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg (-r|-n arg)" << std::endl;
}

template< class ClusterSpectralSequenceT >
void compute_css( SessionConfig conf, int argc, char** argv )
{
    std::ofstream ofs;
    std::string filename = std::string("./results/") + std::string(argv[0]);
    
    for( int i = 1; i < argc; ++i )
    {
        filename += std::string("_") + std::string(argv[i]);
    }
    ofs.open( filename );
    
    // Print status message.
    std::cout << kappa_version( argc, argv )
              << std::endl
              << "------------  Performing computations with the following parameters   ------------" << std::endl
              << "homological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    ofs       << kappa_version( argc, argv )
              << std::endl
              << "------------  Performing computations with the following parameters   ------------" << std::endl
              << "homological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    
    std::cout << "-------- Constructing bases --------" << std::endl;
    ofs       << "-------- Constructing bases --------" << std::endl;
    
    //
    // Compute all bases.
    //
    Clock measure_duration;
    std::cout << "Constructing bases";
    ofs       << "Constructing bases";
    
    ClusterSpectralSequenceT cluster_spectral_sequence( conf.genus, conf.num_punctures, conf.sgn_conv, conf.num_threads, conf.num_remaining_threads );
    typename ClusterSpectralSequenceT::CSSHomologyType homology_E0;
    typename ClusterSpectralSequenceT::CSSHomologyType homology_E1;
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl
              << std::endl;
    ofs       << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl
              << std::endl;
    
    //
    // Print Dimensions of the E^0-page
    //
    std::cout << "Dimensions of the E^0 page" << std::endl;
    ofs       << "Dimensions of the E^0 page" << std::endl;
    for( const auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        const auto& p = basis_it.first;
        const auto& l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        for( int32_t l = 1; l < l_bases.begin()->first ; ++l )
        {
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
        }
        for( const auto& l_basis_it : l_bases )
        {
            const auto& l = l_basis_it.first;
            if( l == 0 )
            {
                continue;
            }
            const auto& cur_basis = l_basis_it.second;
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << cur_basis.size() << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << cur_basis.size() << ";    ";
        }
        std::cout << std::endl;
        ofs       << std::endl;
    }
    
    //
    // Compute E^1 and E^2 terms.
    //
    std::cout << std::endl
              << "-------- Computing E^1-term & E^2-term --------" << std::endl;
    ofs       << std::endl
              << "-------- Computing E^1-term & E^2-term  --------" << std::endl;
    
    // Compute all differentials and homology consecutively.
    for( const auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        const auto& p = basis_it.first;
        const auto& l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        cluster_spectral_sequence.erase_d0();
        cluster_spectral_sequence.erase_d1();
        typename ClusterSpectralSequenceT::MatrixType::DiagonalType diagonal;
        for( const auto& l_basis_it : l_bases )
        {
            const auto& l = l_basis_it.first;
            atomic_uint& current_rank = cluster_spectral_sequence.diff_complex.get_diagonalizer().current_rank;
            uint32_t max_possible_rank(0);
            
            //
            // Construct the first stage of d_1
            //
            std::cout << "Constructing the differential d^1_{" << p << "," << l << "}.";
            ofs       << "Constructing the differential d^1_{" << p << "," << l << "}.";
            measure_duration = Clock();
            cluster_spectral_sequence.gen_d1_stage_1(p,l);
            std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            ofs       << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            
            //
            // Forget the old d0
            //
            cluster_spectral_sequence.erase_d0();
            
            //
            // Generate the new d0
            //
            std::cout << "Constructing the differential d^0_{" << p << "," << l << "}.";
            ofs       << "Constructing the differential d^0_{" << p << "," << l << "}.";
            measure_duration = Clock();
            cluster_spectral_sequence.gen_d0(p,l);
            std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            ofs       << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
            
            //
            // Compute the induced homology of d0.
            //
            measure_duration = Clock(); // Measure duration.
            atomic_uint state(0);       // Set state to 1 iff kernel and torsion are computed. This is done to terminate the 'monitoring thread'.

            // Diagonalzing thread.
            auto partial_homology_thread = std::async( std::launch::async, [&]() -> HomologyField
            {
                auto ret = cluster_spectral_sequence.diff_complex.compute_current_kernel_and_torsion( l );
                state = 1;
                return ret;
            } );

            // Monitoring thread.
            max_possible_rank = std::min( cluster_spectral_sequence.diff_complex.get_current_differential().size1(), cluster_spectral_sequence.diff_complex.get_current_differential().size2() );
            if( (uint32_t)homology_E0[p-1].get_kern(l) > 0 )
            {
                max_possible_rank = std::min( max_possible_rank, (uint32_t)homology_E0[p-1].get_kern(l) );
            }
            auto monitor_thread = std::async( std::launch::async, [&]()
            {
                while( state != 1 )
                {
                    std::cout << "Diagonalization " << current_rank << "/" << max_possible_rank << "\r";
                    std::cout.flush();
                    std::this_thread::sleep_for(std::chrono::milliseconds(50));
                }
            } );

            // Wait for threads to terminate.
            auto partial_homology = partial_homology_thread.get();
            monitor_thread.get();

            //
            // Save results.
            //
            homology_E0[p].set_kern( l, partial_homology.get_kern(l) );
            homology_E0[p-1].set_tors( l, partial_homology.get_tors(l-1) ); // Observe that the computed torsion is located in position l-1.
            
            //
            // Print status message.
            //
            std::cout << "    dim(E^1_{" << (int32_t)(p-1) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E0[p-1].get_kern(l) - homology_E0[p-1].get_tors(l))
                      << "; dim(im d^0_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E0[p-1].get_tors(l))
                      << "; dim(ker d^0_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E0[p].get_kern(l))
                      << "; Duration: " << measure_duration.duration() << " seconds." << std::endl;
            ofs       << "    dim(E^1_{" << (int32_t)(p-1) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E0[p-1].get_kern(l) - homology_E0[p-1].get_tors(l))
                      << "; dim(im d^0_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E0[p-1].get_tors(l))
                      << "; dim(ker d^0_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E0[p].get_kern(l))
                      << "; Duration: " << measure_duration.duration() << " seconds." << std::endl;   
            
            //
            // Save diagonal for the next d1
            //
            diagonal.clear();
            diagonal.swap( cluster_spectral_sequence.diff_complex.get_current_differential().diagonal );
            
            //
            // Detect superflous rows.
            //
            std::list< size_t > ommit_rows;
            for( const auto& it : diagonal )
            {
                ommit_rows.emplace_back( it.first );
            }
            ommit_rows.sort();
            
            //
            // Clean up
            //
            auto& diagonalizer = cluster_spectral_sequence.diff_complex.get_diagonalizer();
            diagonalizer.ommit_rows.clear();
            diagonalizer.ommit_rows.swap(ommit_rows);
            diagonalizer.def = 0;
            diagonalizer.rnk = 0;
            current_rank = 0;
            
            //
            // Compute the induced homology of d1.
            //
            measure_duration = Clock(); // Measure duration.
            state = 0;       // Set state to 1 iff kernel and torsion are computed. This is done to terminate the 'monitoring thread'.
            // Diagonalzing thread.
            cluster_spectral_sequence.prepare_d1_diag();
            partial_homology_thread = std::async( std::launch::async, [&]() -> HomologyField
            {
                auto ret = cluster_spectral_sequence.diff_complex.compute_current_kernel_and_torsion( l );
                state = 1;
                return ret;
            } );

            // Monitoring thread.
            max_possible_rank = std::min( cluster_spectral_sequence.diff_complex.get_current_differential().sec_size1(), cluster_spectral_sequence.diff_complex.get_current_differential().sec_size2() );
            if( (uint32_t)homology_E0[p-1].get_kern(l-1) > 0 )
            {
                max_possible_rank = std::min( max_possible_rank, (uint32_t)(homology_E0[p-1].get_kern(l-1) - homology_E0[p-1].get_tors(l-1) ) );
            }
            monitor_thread = std::async( std::launch::async, [&]()
            {
                while( state != 1 )
                {
                    std::cout << "Diagonalization " << current_rank << "/" << max_possible_rank << "\r";
                    std::cout.flush();
                    std::this_thread::sleep_for(std::chrono::milliseconds(50));
                }
            } );

            // Wait for threads to terminate.
            partial_homology = partial_homology_thread.get();
            monitor_thread.get();
            
            //
            // Save results.
            //
            homology_E1[p].set_kern( l, partial_homology.get_kern(l) );
            homology_E1[p-1].set_tors( l-1, partial_homology.get_tors(l-1) ); // Observe that the computed torsion is located in position l-1.

            //
            // Print status message.
            //
            std::cout << "    dim(E^2_{" << (int32_t)(p-1) << "," << (int32_t)(l-1) << "}) = " << std::setw(4) << (int32_t)(homology_E1[p-1].get_kern(l-1) - homology_E1[p-1].get_tors(l-1) - homology_E0[p-1].get_tors(l-1)) 
                      << "; dim(im d^1_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E1[p-1].get_tors(l-1))
                      << "; dim(ker d^1_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E1[p].get_kern(l))
                      << "; Duration: " << measure_duration.duration() << " seconds." << std::endl;
            ofs       << "    dim(E^2_{" << (int32_t)(p-1) << "," << (int32_t)(l-1) << "}) = " << std::setw(4) << (int32_t)(homology_E1[p-1].get_kern(l-1) - homology_E1[p-1].get_tors(l-1) - homology_E0[p-1].get_tors(l-1)) 
                      << "; dim(im d^1_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E1[p-1].get_tors(l-1))
                      << "; dim(ker d^1_{" << (int32_t)(p) << "," << (int32_t)(l) << "}) = " << std::setw(4) << (int32_t)(homology_E1[p].get_kern(l))
                      << "; Duration: " << measure_duration.duration() << " seconds." << std::endl;
            
            //
            // Forget d1
            //
            cluster_spectral_sequence.erase_d1();
            
            //
            // Restore Diagonal and Diagoanlizer
            //
            diagonal.swap( cluster_spectral_sequence.diff_complex.get_current_differential().diagonal );
            diagonalizer.ommit_rows.clear();
            diagonalizer.def = 0;
            diagonalizer.rnk = 0;
        }
    }

    
    // Print status message.
    std::cout << std::endl
              << "------------  Results   ------------" << std::endl
              << "homological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    ofs       << std::endl
              << "------------  Results   ------------" << std::endl
              << "homological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    
    // Print Dimensions of the E^0-page
    std::cout << "Dimensions of the E^0 page" << std::endl;
    ofs       << "Dimensions of the E^0 page" << std::endl;
    for( const auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        const auto& p = basis_it.first;
        const auto&l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        for( int32_t l = 1; l < l_bases.begin()->first ; ++l )
        {
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
        }
        for( const auto& l_basis_it : l_bases )
        {
            const auto& l = l_basis_it.first;
            if( l == 0 )
            {
                continue;
            }
            const auto& cur_basis = l_basis_it.second;
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << cur_basis.size() << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << cur_basis.size() << ";    ";
        }
        std::cout << std::endl;
        ofs       << std::endl;
    }
    
    // Print Dimensions of the E^1-page
    std::cout << std::endl
              << "Dimensions of the E^1 page" << std::endl;
    ofs       << std::endl
              << "Dimensions of the E^1 page" << std::endl;
    for( const auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        const auto& p = basis_it.first;
        const auto& l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        for( int32_t l = 1; l < l_bases.begin()->first ; ++l )
        {
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
        }
        for( const auto& l_basis_it : l_bases )
        {
            const auto& l = l_basis_it.first;
            if( l == 0 )
            {
                continue;
            }
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << (int32_t)(homology_E0[p].get_kern(l) - homology_E0[p].get_tors(l)) << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << (int32_t)(homology_E0[p].get_kern(l) - homology_E0[p].get_tors(l)) << ";    ";
        }
        std::cout << std::endl;
        ofs       << std::endl;
    }
    
    // Print Dimensions of the E^2-page
    std::cout << std::endl
              << "Dimensions of the E^2 page" << std::endl;
    ofs       << std::endl
              << "Dimensions of the E^2 page" << std::endl;
    for( const auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        const auto& p = basis_it.first;
        const auto& l_bases = basis_it.second.basis;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        for( int32_t l = 1; l < l_bases.begin()->first ; ++l )
        {
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << "0" << ";    ";
        }
        for( const auto& l_basis_it : l_bases )
        {
            const auto& l = l_basis_it.first;
            if( l == 0 )
            {
                continue;
            }
            std::stringstream stream;
            stream << "(" << (int32_t)(p) << "," << (int32_t)(l) << ") =";
            std::cout << std::setw(9) << stream.str() << std::setw(6) << (int32_t)(homology_E1[p].get_kern(l) - homology_E1[p].get_tors(l) - homology_E0[p].get_tors(l)) << ";    ";
            ofs       << std::setw(9) << stream.str() << std::setw(6) << (int32_t)(homology_E1[p].get_kern(l) - homology_E1[p].get_tors(l) - homology_E0[p].get_tors(l)) << ";    ";
        }
        std::cout << std::endl;
        ofs       << std::endl;
    }
    
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
    if( conf.num_threads < 2 )
    {
        conf.num_threads = 2;
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
    
    create_working_directories();
    
    // We may start with the computations.
    if(conf.rational == true)
    {
        compute_css< ClusterSpectralSequenceQ >( conf, argc, argv );
    }
    else if (conf.prime == 2)
    {
        compute_css< ClusterSpectralSequenceBool > ( conf, argc, argv );
    }
    else
    {
        compute_css< ClusterSpectralSequenceZm >( conf, argc, argv );
    }
    
    return 0;
}

