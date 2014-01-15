#include "sessionconfig.hpp"

SessionConfig::SessionConfig(int argc, char **argv) : genus(0), num_punctures(0), rational(0), prime(2), start_p(0), num_threads(0), valid(false)
{
    if( argc >= 4 )
    {
        genus = atoi(argv[1]);
        num_punctures = atoi(argv[2]);
        rational = (atoi(argv[3]) == 0 ? true : false);
        prime = atoi(argv[3]);
        start_p = (argc >= 5 ? atoi(argv[4]) : 0);
        num_threads = (argc >= 6 ? atoi(argv[5]) : 0);
        valid = true;
    }
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
