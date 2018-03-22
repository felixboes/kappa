// The software kappa is a collection of programs to compute the homology of
// the moduli space of surfaces using the radial model.
// Copyright (C) 2013 - 2018  Felix Boes and Anna Hermann
// 
// This file is part of kappa.
// 
// kappa is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// kappa is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with kappa.  If not, see <http://www.gnu.org/licenses/>.


#include <future>
#include <fstream>
#include <iostream>

#include <boost/range/adaptor/reversed.hpp>

#include "kappa.hpp"

void print_usage(int, char** argv)
{
    std::cout << "Usage: " << argv[0] << " -g arg -m arg" << std::endl;
}

void compute_homchain( SessionConfig conf, int argc, char** argv )
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
              << "Saving the Ehrenfried complex associated with the radial model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << std::endl;
    ofs       << kappa_version( argc, argv )
              << std::endl
              << "------------  Performing computations with the following parameters   ------------" << std::endl
              << "Saving the Ehrenfried complex associated with the radial model" << std::endl
              << "genus = " << conf.genus << " punctures = " << conf.num_punctures << std::endl;

    std::cout << "-------- Constructing bases --------" << std::endl;
    ofs       << "-------- Constructing bases --------" << std::endl;

    // Compute all bases.
    Clock measure_duration;
    std::cout << "Constructing bases";
    std::cout.flush();
    ofs << "Constructing bases";

    MonoComplexZStorageOnly monocomplex( conf.genus, conf.num_punctures, conf.sgn_conv, conf.num_threads, conf.num_remaining_threads );
    std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    std::cout << std::endl;
    std::cout.flush();
    ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
    ofs << std::endl;

    std::cout << "-------- Computing Homchain --------" << std::endl;
    ofs << "-------- Computing Homchain --------" << std::endl;

    // Compute all differentials and homology consecutively.
    for( auto& it : boost::adaptors::reverse( monocomplex.basis_complex ) )
    {
        int32_t p = it.first+1;
        if ( p < conf.start_p || p > conf.end_p + 2 )
        {
            continue;
        }

        // Generate a single differential.
        std::cout << "Constructing the " << p << "-th differential of size " << ( monocomplex.basis_complex.count(p) ? monocomplex.basis_complex.at(p).size() : 0 ) << " x " << ( monocomplex.basis_complex.count(p-1) ? monocomplex.basis_complex.at(p-1).size() : 0 );
        std::cout.flush();
        ofs << "Constructing the " << p << "-th differential of size " << ( monocomplex.basis_complex.count(p) ? monocomplex.basis_complex.at(p).size() : 0 ) << " x " << ( monocomplex.basis_complex.count(p-1) ? monocomplex.basis_complex.at(p-1).size() : 0 );

        measure_duration = Clock();
        monocomplex.gen_differential(p);
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;

        measure_duration = Clock();
        std::cout << "Constructing Homchain";
        std::cout.flush();
        ofs << "Constructing Homchain";
        monocomplex.homchain(p);
        int32_t maxdimension = 1 + monocomplex.basis_complex.rbegin()->first;
        monocomplex.homchain(p, true, maxdimension);
        std::cout << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
        ofs << " done. Duration: " << measure_duration.duration() << " seconds." << std::endl;
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

    if( ! ( conf.option_set( "gen" ) && conf.option_set( "pun" ) ) )
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
        SymGrpTuple::parallel_case();
    }
    else
    {
        SymGrpTuple::radial_case();
    }

    create_working_directories();

    compute_homchain( conf, argc, argv );

    return 0;
}
