#ifndef DOUBLECOMPLEXBASIS_HPP
#define DOUBLECOMPLEXBASIS_HPP

#include <iostream>
#include <unordered_set>
#include <string>

#include "libhomology/serialization.hpp"

#include "tuple.hpp"

/**
    The MonoBasis keeps track of the basis elements of a module in a MonoComplex.
**/
struct DoubleComplexBasis
{
    DoubleComplexBasis();
    
    /// Add a basis element.
    uint32_t add_basis_element ( Tuple t );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& mb);
    
    /// Returns the number of basis elements.
    uint64_t size() const;

    /// Returns the index of the Tuple that is stored in the MonoBasis or -1.
    int64_t id_of( const Tuple& t ) const;
    
    /// Stores the orderd basis.
    std::unordered_set< Tuple, HashTuple > basis;
    
    /// Returns a const reference to the container.
    const std::unordered_set< Tuple, HashTuple >& get_container() const;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a MonoBasis to file.
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
    
    /// boost::serialization method that we use to load a MonoBasis from file.
    template<class Archive>
    void load(Archive & ar, const unsigned int)
    {
        size_t size;
        Tuple t;
        
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

std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& basis);

#endif // DOUBLECOMPLEXBASIS_HPP
