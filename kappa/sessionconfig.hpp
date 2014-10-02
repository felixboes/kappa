#ifndef SESSIONCONFIG_HPP
#define SESSIONCONFIG_HPP

#include <limits>
#include <string>

#include <boost/program_options.hpp>

#include "libhomology/homology.hpp"

/**
 *  \brief The sign convention in the Ehrenfried complex.
 *  Here, all_signs corresponds to the orientation coefficients,
 *  no_orientation_signs corresponds to the constant coefficient system and
 *  no_signs is allowed to be used iff the coeffiencts are a \f$ \mathbb F_2 \f$-algebra.
 */ 
enum SignConvention
{
    all_signs,              ///< use the orientation coefficients
    no_orientation_sign,    ///< use the constant coefficients system
    no_signs                ///< use the constant coefficients system and assume that we use a \f$ \mathbb F_2\f$ algebra
};

/**
 *  The SessionConfig struct is used derive a configuration from a given set of command line arguments.
 *  Doing so, we can write several main()s more compact.
 */
struct SessionConfig
{
    SessionConfig() :
        genus(0),
        num_punctures(0),
        rational(0),
        prime(2),
        parallel(false),
        num_threads(1),
        num_remaining_threads(0),
        start_p(0),
        end_p(std::numeric_limits<uint32_t>::max()),
        sgn_conv(all_signs),
        create_cache(true),
        valid(false),
        print_help(false)
    {}
    
    SessionConfig( const int argc, char** argv );
    
    /// Sometimes we have to work a litte, before we can use a configurration. Here we check parameters and setup Z_p coefficients.
    bool setup_configuration();
    
    /// @returns true iff the command line option was set.
    bool option_set( const std::string opt ) const;
    
    /// @returns wether the given configuration is valid or not.
    operator bool() const {return valid;}
    
    boost::program_options::options_description desc;
    boost::program_options::variables_map vm;
    
    uint32_t genus;     ///< The genus of the surface.
    uint32_t num_punctures; ///< The number of punctures of the surface.
    bool rational;      ///< True iff we work over the field \f$ \mathbb Q\f$.
    uint32_t prime;     ///< The prime \f$ s \f$ if we work over $\mathbb Z / s \mathbb Z$.
    bool parallel;      ///< True iff we use the parallel model.
    uint32_t num_threads;   ///< The number of total threads in use.
    uint32_t num_remaining_threads; ///< The number of threads that collect work in the diagonalizing process.
    int32_t start_p;        ///< The homological degree we want to start at.
    int32_t end_p;          ///< The homological degree we want to end with.
    int32_t first_basis;    ///< The first basis we want to consider.
    int32_t last_basis;     ///< The last basis we want to consider.
    SignConvention sgn_conv;    ///< The sign convention we use.
    bool create_cache;      ///< Stores diagonalized matrices iff this is true.
    bool valid;             ///< Is set to true if a given configuration is valid.
    bool print_help;        ///< Prints the help to screen iff this is true.
};

#endif // SESSIONCONFIG_HPP
