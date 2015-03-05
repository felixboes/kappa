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
    void add_basis_element ( HighCell t );
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& mb);
    
    /// Returns the number of basis elements.
    uint64_t size_red() const;
    uint64_t size_col() const;
    uint64_t size_ess() const;

    /// Generates the enumeration of the basis elements.
    /// This should be called as soon as the basis is completely generated.
    void generate_indices();

    /// Returns the index of the HighCell that is stored in the Basis or -1.
    int64_t id_of( const HighCell& t ) const;
    
    /// Returns the index of the HighCell associated to the Tuple that is stored in the Basis or -1.
    int64_t id_of( const Tuple& t ) const;
    
    /// Stores the orderd basis.
    std::set< HighCell > basis_red;
    std::set< HighCell > basis_col;
    std::set< HighCell > basis_ess;
    
    /// Returns a const reference to the container.
    const std::set< HighCell >& get_container_red() const;
    const std::set< HighCell >& get_container_col() const;
    const std::set< HighCell >& get_container_ess() const;
    
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int) const
    {
        ar & basis_red & basis_col & basis_ess;
    }
};

std::ostream& operator<< (std::ostream& stream, const DoubleComplexBasis& basis);

#endif // DOUBLECOMPLEXBASIS_HPP
