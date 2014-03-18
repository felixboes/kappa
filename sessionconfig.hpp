#ifndef SESSIONCONFIG_HPP
#define SESSIONCONFIG_HPP

#include <limits>

#include <homology.hpp>

    
enum SignConvention
{
    all_signs,
    no_orientation_sign,
    no_signs
};

/**
 *  The SessionConfig struct is used derive a configuration from a given set of command line arguments.
 *  Doing so, we can write several main()s more compact.
 */
struct SessionConfig
{
    SessionConfig() : genus(0), num_punctures(0), rational(0), prime(2), start_p(0), end_p(std::numeric_limits<uint32_t>::max()), sgn_conv(all_signs), valid(false) {}
    SessionConfig( int argc, char** argv );
    
    /// Sometimes we have to work a litte, before we can use a configurration. Here we check parameters and setup Z_p coefficients.
    bool setup_configuration();
    
    /// @returns wether the given configuration is valid or not.
    inline operator bool() {return valid;}
    
    uint32_t genus;
    uint32_t num_punctures;
    bool rational;
    uint32_t prime;
    uint32_t num_threads;
    uint32_t start_p;
    uint32_t end_p;
    SignConvention sgn_conv;
    bool valid;
};

#endif // SESSIONCONFIG_HPP
