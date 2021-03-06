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


#ifndef SESSIONCONFIG_HPP
#define SESSIONCONFIG_HPP

#include <limits>
#include <string>

#include <boost/program_options.hpp>

#include <libhomology/homology.hpp>

#include "misc.hpp"

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
    SessionConfig();
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
    bool apply_base_changes;    ///< Apply base changes.
    bool create_cache;      ///< Stores diagonalized matrices iff this is true.
    int32_t limit_mem;      ///< Stores the percentage of the total memory available that is allowed to be used.
    bool valid;             ///< Is set to true if a given configuration is valid.
    bool print_help;        ///< Prints the help to screen iff this is true.
};

#endif // SESSIONCONFIG_HPP
