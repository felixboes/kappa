#ifndef DOUBLECOMPLEXBASIS_HPP
#define DOUBLECOMPLEXBASIS_HPP

#include <iostream>
#include <unordered_set>
#include <string>

#include "libhomology/serialization.hpp"

#include "high_cell.hpp"

/**
    The MonoBasis keeps track of the basis elements of a module in a MonoComplex.
**/
struct DoubleComplexBasis
{
    DoubleComplexBasis();
    
    /// Add a basis element.
    uint32_t add_basis_element ( HighCell t );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& mb);
    
    /// Returns the number of basis elements.
    uint64_t size_h() const;
    uint64_t size_h_1() const;

    /// Returns the index of the HighCell that is stored in the MonoBasis or -1.
    int64_t id_of( const HighCell& t ) const;
    
    /// Stores the orderd basis.
    std::unordered_set< HighCell, HashHighCell > basis_h;
    std::unordered_set< HighCell, HashHighCell > basis_h_1;
    
    /// Returns a const reference to the container.
    const std::unordered_set< HighCell, HashHighCell >& get_container_h() const;
    const std::unordered_set< HighCell, HashHighCell >& get_container_h_1() const;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a MonoBasis to file.
    template<class Archive>
    void save(Archive & ar, const unsigned int) const
    {
        // In order to load an unorderd_set we need to know the exact number of elemets that are stored.
        size_t size = basis_h.size();
        ar & size;
        for( const auto& it : basis_h )
        {
            ar & it;
        }
        size = basis_h_1.size();
        ar & size;
        for( const auto& it : basis_h_1 )
        {
            ar & it;
        }
    }
    
    /// boost::serialization method that we use to load a MonoBasis from file.
    template<class Archive>
    void load(Archive & ar, const unsigned int)
    {
        size_t size;
        HighCell t;
        
        ar & size;
        for( size_t i = 0; i < size; ++i )
        {
            ar & t;
            basis_h.insert(t);
        }
        ar & size;
        for( size_t i = 0; i < size; ++i )
        {
            ar & t;
            basis_h_1.insert(t);
        }
    }
    
    // This is required as saving and loading are different methods.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& basis);

#endif // DOUBLECOMPLEXBASIS_HPP
