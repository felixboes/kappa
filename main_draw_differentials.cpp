#include <iostream>

#include <Magick++.h>

#include "kappa.hpp"

void print_usage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " genus num_punctures (rational:0 | Z_r: r > 0) optional:number of threads for paralellization optional: first and last p for which to compute H_p." << std::endl;
}

template< class ClusterSpectralSequenceT >
void draw_differentials( SessionConfig conf, int argc, char** argv )
{
    std::cout << "-------- Constructing bases --------" << std::endl;
    
    // Compute all bases.
    Clock measure_duration;
    std::cout << "Constructing bases";
    std::cout.flush();
    
    ClusterSpectralSequenceT cluster_spectral_sequence( conf.genus, conf.num_punctures, conf.sgn_conv );
    
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    std::cout << std::endl;
    std::cout.flush();
    
    for( auto& basis_it : cluster_spectral_sequence.basis_complex )
    {
        auto p = basis_it.first;

        if ( p < conf.start_p || p > conf.end_p + 1 )
        {
            continue;
        }
        
        // Generate a single differential.
        std::cout << "Constructing and drawing the differential " << p << "-th differential ";
        std::cout.flush();

        measure_duration = Clock();
        cluster_spectral_sequence.draw_differential(p);
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        std::cout.flush();
    }
}

int main( int argc, char** argv )
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
    
    // Initialize the image magick library
    try {
        Magick::InitializeMagick(argv[0]);
    }
    catch( Magick::Exception &error_ )
    {
        std::cout << "Caught exception: " << error_.what() << std::endl;
        return 1;
    }
    
    
    // We may start with the computations.
    if(conf.rational == true)
    {
        draw_differentials< ClusterSpectralSequenceQ >( conf, argc, argv );
    }
    else
    {
        draw_differentials< ClusterSpectralSequenceZm >( conf, argc, argv );
    }
    
    return 0;
}
