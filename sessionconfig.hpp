#ifndef SESSIONCONFIG_HPP
#define SESSIONCONFIG_HPP

#include <homology.hpp>

struct SessionConfig
{
    SessionConfig() : genus(0), num_punctures(0), rational(0), prime(2), valid(false) {}
    SessionConfig( int argc, char** argv );
    
    bool setup_configuration();
    operator bool() {return valid;}
    
    uint32_t genus;
    uint32_t num_punctures;
    bool rational;
    uint32_t prime;
    bool valid;
};

#endif // SESSIONCONFIG_HPP
