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


#include "sessionconfig.hpp"

SessionConfig::SessionConfig() :
    genus(0),
    num_punctures(0),
    rational(false),
    prime(2),
    parallel(false),
    num_threads(1),
    num_remaining_threads(0),
    start_p(0),
    end_p(0),
    first_basis(0),
    last_basis(0),
    sgn_conv(all_signs),
    apply_base_changes(false),
    create_cache(false),
    limit_mem(0),
    valid(false),
    print_help(false)
{}

SessionConfig::SessionConfig( const int argc, char **argv ) : 
    desc("Command line options"),
    genus(0),
    num_punctures(0),
    rational(false),
    prime(2),
    parallel(false),
    num_threads(1),
    num_remaining_threads(0),
    start_p(0),
    end_p(0),
    first_basis(0),
    last_basis(0),
    apply_base_changes(false),
    create_cache(false),
    limit_mem(0),
    valid(false),
    print_help(false)
{
    
    desc.add_options()
            ("help,h", "produce help message")
            ("gen,g", boost::program_options::value(&genus), "the genus of the Riemann surfaces")
            ("pun,m", boost::program_options::value(&num_punctures), "the number of punctures of the Riemann surfaces")
            ("rat,q", "uses rational numbers as coefficients")
            ("fin,s", boost::program_options::value(&prime), "uses the finite field F_s with s a prime number")
            ("parallel", "uses radial or parallel slit configurations")
            ("num_working_threads,t", boost::program_options::value(&num_threads), "the number of threads used for work in matrix computations")
            ("num_remaining_threads", boost::program_options::value(&num_remaining_threads), "the number of additional threads used in diagonalization")
            ("apply_base_changes", "applies base changes to the matrices.")
            ("cache", "cache the diagonalized matrices")
            ("limit", boost::program_options::value(&limit_mem), "limit used memory to arg% of the total memory available")
            ("first_diff", boost::program_options::value(&start_p), "start with the differential first_diff")
            ("last_diff", boost::program_options::value(&end_p), "end with the differential last_diff")
            ("first_basis", boost::program_options::value(&first_basis), "start with the basis first_basis")
            ("last_basis", boost::program_options::value(&last_basis), "end with the basis last_basis")
    ;

    // Parse command line options
    try{
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        try{
            boost::program_options::notify(vm);
        }
        catch(boost::program_options::required_option& e)
        {
            std::cerr << "Error: " << e.what() << std::endl << std::endl;
            return;
        }
    }
    catch(boost::program_options::error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl << std::endl;
        return;
    }

    // Configure Session
    if( vm.count("limit") == false )
    {
        limit_mem = 0;
    }
    else
    {
        if(limit_mem < 0)
        {
            limit_mem = 0;
        }
        else if( limit_mem > 100 )
        {
            limit_mem = 100;
        }
    }
    if( vm.count("help") == true )
    {
        print_help = true;
    }
    if( vm.count("fin") == false )
    {
        rational = true;
    }
    if( vm.count("last_diff") == false )
    {
        end_p = 4*genus+2*num_punctures; // this is 2h.
    }
    if( vm.count("last_basis") == false )
    {
        last_basis = 4*genus+2*num_punctures; // this is 2h.
    }
    parallel = vm.count("parallel");
    apply_base_changes = vm.count("apply_base_changes");
    create_cache = vm.count("cache");
    
    valid = true;
}

bool SessionConfig::option_set( const std::string opt ) const
{
    // Check if option is set and if was really set.
    return vm.count(opt) && !vm.at(opt).defaulted();
}

bool SessionConfig::setup_configuration()
{
    if( valid )
    {
        if( limit_mem > 0 )
        {
            limit_memory(limit_mem);
        }
        if( num_punctures >= 2 )
        {
            sgn_conv = all_signs;
        }
        else
        {
            sgn_conv = no_orientation_sign;
        }
        
        if( rational == false )
        {
            bool is_prime = true;
            for( uint32_t i = 2; i < prime; ++i )
            {
                if( prime % i == 0 )
                {
                    is_prime = false;
                    break;
                }
            }
            
            if( is_prime )
            {
                Zm::set_modulus(prime);
                if(prime == 2)
                {
                    sgn_conv = no_signs;
                }
            }
            else
            {
                std::cout << "Die Zahl p ist nicht prim." << std::endl;
                return false;
            }
        }
        
        return true;
    }
    else
    {
        return false;
    }
}
