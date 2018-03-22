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


#ifndef CSSBASIS_HPP
#define CSSBASIS_HPP

#include <iostream>
#include <unordered_set>

#include "libhomology/serialization.hpp"

#include "sym_grp_tuple.hpp"

/**
    The CSSBasis keeps track of the basis elements of a module in the cluster spectral sequence.
**/
struct CSSBasis
{
    /// Stores the orderd basis where every cell has the same number l of cluster.
    typedef std::unordered_set< SymGrpTuple, HashSymGrpTuple > LBasisType;
    
    /// Stores the orderd basis, sorted by cluster sizes.
    typedef std::map< int32_t , LBasisType > BasisType;
    
    /// Add a basis element.
    int32_t add_basis_element ( SymGrpTuple t );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const CSSBasis& cssb);
    
    /// Returns the number of basis elements that have cluster size l.
    int32_t size( const int32_t l ) const;
    
    /// Returns the number of basis elements.
    int32_t total_size() const;

    /// Returns the relative (i.e. the cluster-) index of the Tuple that is stored in the CSSBasis or -1.
    int32_t id_of( SymGrpTuple& t ) const;
    
    /// Returns the total index of the Tuple that is stored in the CSSBasis or -1.
    int32_t total_id_of( SymGrpTuple& t ) const;
    
    BasisType basis;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a CSSBasis to file.
    template<class Archive>
    void save(Archive & ar, const unsigned int) const
    {
        // In order to load an unorderd_set we need to know the exact number of elemets that are stored.
        size_t num_cluster_sizes = basis.size(); // this is the number of cluster sizes that occur.
        ar & num_cluster_sizes;
        
        for( const auto& it : basis )
        {
            ar & it.first; // This is a specific l.
            for( const auto& b : it.second )
            {
                ar & b; // This is a basis element with exactly l clusters.
            }
        }
    }
    
    /// boost::serialization method that we use to load a CSSBasis from file.
    template<class Archive>
    void load(Archive & ar, const unsigned int)
    {
        size_t num_cluster_sizes;
        
        ar & num_cluster_sizes;
        for( size_t i = 0; i < num_cluster_sizes; ++i )
        {
            size_t l;
            ar & l;
            SymGrpTuple t;
            for( size_t j = 0; j < l; ++l )
            {
                ar & t;
                basis[l].insert(t);
            }
        }
    }
    
    // This is required as saving and loading are different methods.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

std::ostream& operator<< (std::ostream& stream, const CSSBasis& basis);

#endif // CSSBASIS_HPP
