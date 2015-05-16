#include <iostream>

#include <Magick++.h>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg (-r|-n arg)" << std::endl;
}

template< class MatrixType >
void draw_differential( const MatrixType& matrix, std::string filename )
{
    if( matrix.size1() == 0 || matrix.size2() == 0 )
    {
        return;
    }

    const MatrixType m(1,1);
    const auto zero = m.at(0,0);

    Magick::Image picture(  Magick::Geometry( matrix.size1(), matrix.size2() ), Magick::Color("white") );
    try{
        // Define color and thickness of a point.
        picture.strokeColor("red");
        picture.fillColor("red");

        // Draw Column
        for( size_t j = 0; j < matrix.size2(); ++j )
        {
            for( size_t i = 0; i < matrix.size1(); ++i )
            {
                if( matrix.at(i,j) != zero )
                {
                    picture.draw( Magick::DrawableCircle(i, j, i+1, j ) );
                }
            }
        }
        picture.write( filename );
    }
    catch( Magick::Exception &error_ )
    {
        std::cout << "Caught exception: " << error_.what() << std::endl;
    }
}

template< class MonoComplexT >
void draw_differentials( SessionConfig conf, int argc, char** argv )
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
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    std::cout << std::endl;
    std::cout.flush();
    ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    ofs << std::endl;

    std::cout << "-------- Drawing Differentials --------" << std::endl;
    ofs << "-------- Drawing Differentials --------" << std::endl;

    // Compute all differentials and homology consecutively.
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
        monocomplex.gen_differential(p+1);
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;

        std::string filename = "results/Differential_";
        filename += std::to_string(conf.genus) + std::string("_") + std::to_string(conf.num_punctures) + std::string("_-_") + std::to_string(p) + std::string(".png");
        draw_differential( monocomplex.get_current_differential(), filename );
    }
}

int main( int argc, char** argv )
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
    
    if ( conf.setup_configuration() == false )
    {
        std::cout << "The configuration could not been setup." << std::endl;
        return 2;
    }
    
    // Initialize the image magick library
    try {
        Magick::InitializeMagick(argv[0]);
    }
    catch( Magick::Exception &error_ )
    {
        std::cout << "Caught exception: " << error_.what() << std::endl;
        return 1;
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

    if(conf.rational == true)
    {
        draw_differentials< MonoComplexQ >( conf, argc, argv );
    }
    else
    {
        draw_differentials< MonoComplexZm >( conf, argc, argv );
    }
    
    return 0;
}
