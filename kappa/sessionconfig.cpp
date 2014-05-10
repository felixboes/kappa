#include "sessionconfig.hpp"

SessionConfig::SessionConfig( int argc, char **argv ) : 
    desc("Command line options"),
    genus(0),
    num_punctures(0),
    rational(false),
    prime(2),
    num_threads(1),
    start_p(0),
    end_p(0),
    first_basis(0),
    last_basis(0),
    valid(false),
    print_help(false)
{
    
    desc.add_options()
            ("help,h", "produce help message")
            ("gen,g", boost::program_options::value(&genus), "the genus of the Riemannian surfaces")
            ("pun,m", boost::program_options::value(&num_punctures), "the number of punctures of the Riemannian surfaces")
            ("rat,r", "uses rational numbers as coefficients")
            ("fin,n", boost::program_options::value(&prime), "uses the finite field F_n with n a prime number")
            ("num_threads,t", boost::program_options::value(&num_threads), "the number of threads used in matrix computations")
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
            Zm::set_modulus(prime);
            if(prime == 2)
            {
                sgn_conv = no_signs;
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}