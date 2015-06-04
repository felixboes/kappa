#include <future>
#include <fstream>
#include <iostream>

#include <boost/range/adaptor/reversed.hpp>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg (-q|-n arg)" << std::endl;
}

template< class MonoComplexT >
void compute_homology( SessionConfig conf, int argc, char** argv )
{
    // Prepare status messages.
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
              << "Cohomological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    ofs       << kappa_version( argc, argv )
              << std::endl
              << "------------  Performing computations with the following parameters   ------------" << std::endl
              << "Cohomological Ehrenfried complex associated with the " << (conf.parallel == true ? "parallel" : "radial") << " model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << " coefficients = " << ( conf.rational == true ? "Q" : ("Z/" + std::to_string(conf.prime) + "Z") ) << std::endl
              << std::endl;
    
    std::cout << "-------- Constructing bases --------" << std::endl;
    ofs       << "-------- Constructing bases --------" << std::endl;
    
    // Compute all bases.
    Clock measure_duration;
    std::cout << "Constructing bases";
    std::cout.flush();
    ofs << "Constructing bases";

    MonoComplexT monocomplex( conf.genus, conf.num_punctures, conf.sgn_conv, conf.num_threads, conf.num_remaining_threads );
    typename MonoComplexT::HomologyType homology;
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    std::cout << std::endl;
    std::cout.flush();
    ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    ofs << std::endl;

    std::cout << "-------- Computing Homology --------" << std::endl;
    ofs << "-------- Computing Homology --------" << std::endl;

    // Compute all differentials and homology consecutively.
    SagemathInterface sage;
    for( auto& it : boost::adaptors::reverse( monocomplex.bases ) )
    {
        int32_t p = it.first;
        if ( p < conf.start_p || p > conf.end_p + 2 )
        {
            continue;
        }

        // Generate a single differential.
        std::cout << "Constructing the " << p << "-th differential of size " << monocomplex.bases.at(p).size() << " x " << ( monocomplex.bases.count(p+1) ? monocomplex.bases.at(p+1).size() : 0 );
        std::cout.flush();
        ofs << "Constructing the " << p << "-th differential of size " << monocomplex.bases.at(p).size() << " x " << ( monocomplex.bases.count(p+1) ? monocomplex.bases.at(p+1).size() : 0 );

        measure_duration = Clock();
        monocomplex.gen_differential_sage(p+1, sage);
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;

        // Diagoanlize differetnial and save results.
        measure_duration = Clock();
        sage.compute_rank();

        // Print status message.
        std::cout << "Diagonalization done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
        ofs << "Diagonalization done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << std::endl;
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
//    else if (conf.prime == 2)
//    {
//        compute_homology< MonoComplexBool > ( conf, argc, argv);
//    }
//    else
//    {
//        compute_homology< MonoComplexZm >( conf, argc, argv );
//    }

    return 0;
}
