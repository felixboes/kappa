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
    std::string filename = std::string("./results/") + std::string(argv[0]);

    for( int i = 1; i < argc; ++i )
    {
        filename += std::string("_") + std::string(argv[i]);
    }

    ofs.open( filename );
    
    // Print status message.
    std::cout << std::endl
              << "Program name and parameters: ";
    for( int i = 0; i < argc; ++i )
    {
        std::cout << argv[i];
        if( i+1 < argc )
        {
            std::cout << " ";
        }
    }         
    std::cout << std::endl
              << "Program version: " << program_version_by_git << std::endl
              << "Date: " << current_date() << std::endl
              << std::endl
              << "------------  Performing computations with the following parameters   ------------" << std::endl
              << "homological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    ofs       << std::endl
              << "Program name and parameters: ";
    for( int i = 0; i < argc; ++i )
    {
        ofs << argv[i];
        if( i+1 < argc )
        {
            ofs << " ";
        }
    }
    ofs       << std::endl
              << "Program version: " << program_version_by_git << std::endl
              << "Date: " << current_date() << std::endl
              << std::endl
              << "------------  Performing computations with the following parameters   ------------" << std::endl
              << "homological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    
    std::cout << "-------- Constructing bases --------" << std::endl;
    ofs       << "-------- Constructing bases --------" << std::endl;

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
    
    if( conf.create_cache == true ) // Save all bases elements
    {
        std::string path_prefix = 
            std::string("./cache/bases_") + std::string( (conf.parallel == true ? "parallel" : "radial") ) + std::string("/") + 
            std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_";
        std::string check_writable_prefix = 
            std::string("./cache/list_of_files_that_should_not_be_overwritten/bases_") +
            std::string( (conf.parallel == true) ? "parallel_" : "radial_" ) + std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_";
        
        for( const auto& it : monocomplex.bases )
        {
            const auto& p = it.first;
            
            // Check if file may be overwritten
            if( file_exists(check_writable_prefix + std::to_string(p) ) == true )
            {
                std::cout << "Do not overwrite " << path_prefix + std::to_string(p) << std::endl;
            }
            else
            {
                // Store the p-th basis.
                save_to_file_bz2<MonoBasis>(it.second, path_prefix + std::to_string(p));
                touch( check_writable_prefix + std::to_string(p) );
            }
        }
    }

    std::cout << "-------- Computing Homology --------" << std::endl;
    ofs << "-------- Computing Homology --------" << std::endl;

    // Compute all differentials and homology consecutively.
    for( auto& it : monocomplex.bases )
    {
        int32_t p = it.first;
        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }

        // Generate a single differential.
        std::cout << "Constructing the " << p << "-th differential of size " << monocomplex.bases[p].size() << " x " << monocomplex.bases[p-1].size();
        std::cout.flush();
        ofs << "Constructing the " << p << "-th differential of size " << monocomplex.bases[p].size() << " x " << monocomplex.bases[p-1].size();

        measure_duration = Clock();
        monocomplex.gen_differential(p);

        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;

        // Diagoanlize differetnial and save results.
        measure_duration = Clock();
        uint32_t max_possible_rank( std::min( monocomplex.num_rows(), monocomplex.num_cols() ) );
        if( (uint32_t)homology.get_kern(p-1) > 0 )
        {
            max_possible_rank = std::min( max_possible_rank, (uint32_t)homology.get_kern(p-1) );
        }
        auto partial_homology = monocomplex.diagonalize_current_differential( p, max_possible_rank, false );
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
        
        if( conf.create_cache == true ) // Save base changes and triangular shape.
        {
            // The prefix is cache/differentials_(parallel|radial)/(coeff)_(genus)_(punctures)_(p)
            const auto & cur_differential = monocomplex.get_current_differential();
            std::string path_prefix =
                "./cache/differentials_" + std::string( ( conf.parallel == true ? "parallel" : "radial") ) + std::string("/") +
                (conf.rational == true ? std::string("q") : std::string("s") + std::to_string(conf.prime) ) + std::string("_") +
                std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_" + std::to_string(p) + std::string("_");
            std::string check_writable_prefix =
                "./cache/list_of_files_that_should_not_be_overwritten/differentials_" + std::string( ( conf.parallel == true ? "parallel" : "radial") ) + std::string("_") +
                (conf.rational == true ? std::string("q") : std::string("s") + std::to_string(conf.prime) ) + std::string("_") +
                std::to_string(conf.genus) + "_" + std::to_string(conf.num_punctures) + "_" + std::to_string(p) + std::string("_");
            
            if( file_exists(check_writable_prefix + "diagonal" ) == true )
            {
                std::cout << "Do not overwrite " << path_prefix << "diagonal" << std::endl;
            }
            else
            {
                cur_differential.cache_diagonal( path_prefix + "diagonal" );
                touch( check_writable_prefix + "diagonal" );
            }
            
            if( file_exists(check_writable_prefix + "base_changes" ) == true )
            {
                std::cout << "Do not overwrite " << path_prefix << "base_changes" << std::endl;
            }
            else
            {
                cur_differential.cache_base_changes( path_prefix + "base_changes" );
                touch( check_writable_prefix + "base_changes" );
            }
            
            if( file_exists(check_writable_prefix + "triangular" ) == true )
            {
                std::cout << "Do not overwrite " << path_prefix << "triangular" << std::endl;
            }
            else
            {
                cur_differential.cache_triangular_shape( path_prefix + "triangular" );
                touch( check_writable_prefix + "triangular" );
            }
        }
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
    
    create_working_directories();

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
