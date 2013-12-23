#include "sessionconfig.hpp"

SessionConfig::SessionConfig(int argc, char **argv) : genus(0), num_punctures(0), rational(0), prime(2), valid(false)
{
    if( argc >= 4)
    {
        genus = atoi(argv[1]);
        num_punctures = atoi(argv[2]);
        rational = (atoi(argv[3]) == 0 ? true : false);
        prime = atoi(argv[3]);
        valid = true;
    }
}

bool SessionConfig::setup_configuration()
{
    if( valid )
    {
        if( rational == false )
        {
            Zm::set_modulus(prime, 1);
        }
        return true;
    }
    else
    {
        return false;
    }
}
