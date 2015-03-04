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
    uint64_t size_red() const;
    uint64_t size_col() const;
    uint64_t size_ess() const;

    /// Returns the index of the HighCell that is stored in the MonoBasis or -1.
    int64_t id_of( const HighCell& t ) const;
    
    /// Stores the orderd basis.
    std::unordered_set< HighCell, HashHighCell > basis_red;
    std::unordered_set< HighCell, HashHighCell > basis_col;
    std::unordered_set< HighCell, HashHighCell > basis_ess;
    
    /// Returns a const reference to the container.
    const std::unordered_set< HighCell, HashHighCell >& get_container_red() const;
    const std::unordered_set< HighCell, HashHighCell >& get_container_col() const;
    const std::unordered_set< HighCell, HashHighCell >& get_container_ess() const;
    
    friend class boost::serialization::access;
    
    /// @warning The serialization library from boost does not yet support unorderd maps (we use boost in the version 1.49). Therefore we must provide a workaround.
    /// boost::serialization method that we use to save a MonoBasis to file.
    template<class Archive>
    void save(Archive & ar, const unsigned int) const
    {
        // In order to load an unorderd_set we need to know the exact number of elemets that are stored.
        size_t size = basis_red.size();
        ar & size;
        for( const auto& it : basis_red )
        {
            ar & it;
        }
        
        size = basis_col.size();
        ar & size;
        for( const auto& it : basis_col )
        {
            ar & it;
        }
        
        size = basis_ess.size();
        ar & size;
        for( const auto& it : basis_ess )
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
            basis_red.insert(t);
        }
        
        ar & size;
        for( size_t i = 0; i < size; ++i )
        {
            ar & t;
            basis_col.insert(t);
        }
        
        ar & size;
        for( size_t i = 0; i < size; ++i )
        {
            ar & t;
            basis_ess.insert(t);
        }
    }
    
    // This is required as saving and loading are different methods.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& basis);

#endif // DOUBLECOMPLEXBASIS_HPP
