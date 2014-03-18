#include "sessionconfig.hpp"

SessionConfig::SessionConfig(int argc, char **argv) : genus(0), num_punctures(0), rational(0), prime(2), num_threads(0), start_p(0), end_p(std::numeric_limits<uint32_t>::max()), valid(false)
{
    if( argc >= 4 )
    {
        genus = atoi(argv[1]);
        num_punctures = atoi(argv[2]);
        rational = (atoi(argv[3]) == 0 ? true : false);
        prime = atoi(argv[3]);
        num_threads = (argc >= 5 ? atoi(argv[4]) : 0);
        start_p = (argc >= 6 ? atoi(argv[5]) : 0);
        end_p = (argc >= 7 ? atoi(argv[6]) : std::numeric_limits<uint32_t>::max() - 100);
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
