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


#ifndef EHR_BASIS_HPP
#define EHR_BASIS_HPP

#include <iostream>
#include <unordered_set>
#include <string>

#include "libhomology/serialization.hpp"
#include "sym_grp_tuple.hpp"

/**
    The EhrBasis keeps track of the basis elements of a module in a EhrComplex.
**/
template <class TupleT>
struct EhrBasis
{
    typedef TupleT TupleType;
    typedef EhrBasis< TupleT > ThisType;

    EhrBasis();
    
    /// Add a basis element.
    uint32_t add_basis_element ( TupleT t );

    /// Add a basis element that is not a multiple of a.
    uint32_t add_basis_element_reduced ( TupleT t );
    
    /// output stream
    template <class T>
    friend std::ostream& operator<< (std::ostream& stream, const EhrBasis< T >& mb);

    /// Returns the number of basis elements.
    uint64_t size() const;

    /// Returns the index of the Tuple that is stored in the EhrBasis or -1.
    int64_t id_of( const TupleT& t ) const;
    
    /// Stores the orderd basis.
    std::unordered_set< TupleT, HashSymGrpTuple > basis;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a EhrBasis to file.
    template<class Archive>
    void save(Archive & ar, const unsigned int) const
    {
        // In order to load an unorderd_set we need to know the exact number of elemets that are stored.
        size_t size = basis.size();
        ar & size;
        for( const auto& it : basis )
        {
            ar & it;
        }
    }
    
    /// boost::serialization method that we use to load a EhrBasis from file.
    template<class Archive>
    void load(Archive & ar, const unsigned int)
    {
        size_t size;
        TupleT t;
        
        ar & size;
        for( size_t i = 0; i < size; ++i )
        {
            ar & t;
            basis.insert(t);
        }
    }
    
    // This is required as saving and loading are different methods.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

template< class TupleT >
EhrBasis< TupleT > load_ehr_basis(const uint32_t g, const uint32_t m, const int32_t p, const bool radial);

template< class TupleT >
std::ostream& operator<< (std::ostream& stream, const EhrBasis< TupleT >& basis);

#endif // EHR_BASIS_HPP
