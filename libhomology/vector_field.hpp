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


#ifndef VECTOR_FIELD_HPP
#define VECTOR_FIELD_HPP

// Description:
//
// This header defines a vector type with coefficients in an (arbitrary) field.
// Such a field must posess operators like 'operator/= ()'.
//
// Moreover this header offers vectors with Q and Zm coefficients, namely VectorQ and VectorZm.

#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "field_coefficients.hpp"
#include "matrix_field.hpp"

/**
 *  This template class defines a vector type subject to the field coeffiecients 'CoefficientT'.
 *  Our implementation mimes the functionality of ublas::matrix but performs faster.
 *  You can access elements via operator() e.g. VectorField v(10) M(0) = CoefficientT(-3);
 */
template < class CoefficientT >
class VectorField
{
public:
    typedef CoefficientT CoefficientType;
    typedef VectorField< CoefficientType > ThisType;
    typedef std::vector< CoefficientT > VectorStorageType;  ///< This realizes the implementation of the data.
    
    /**
     *  Creates a vector with zero entries.
     */
    VectorField();
    
    /**
     *  Creates a vector of a given dimension.
     *  The entries are determined by the standard constructor of CoefficientT.
     *  @warning If the standard constructor of CoefficientT does not create a coefficient with value zero (e.g. usind standard int types)
     *  the result differ from your expectation.
     *  In order to get a zero vector you may use clear();
     */
    VectorField( const size_t dimension );
    
    /**
     *  In order to access elements of the vector $\f v\f$ you want to use this function.
     *  @return The function returns a reference to the given entry.
     *  @todo throw an exception if necessary i.e. if \f$ v_i \f$ is not a valid entry.
     */ 
    CoefficientT & operator()( const size_t i );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    const CoefficientT& at( const size_t i ) const;
    
    /**
     *  adds a given vector to this vector.
     */
    ThisType& operator+=( const ThisType& argument );
    
    /**
     *  substracts a given vector from this vector.
     */
    ThisType& operator-=( const ThisType& argument );
    
    /**
     *  multiply this vector with a given scalar.
     */
    ThisType& operator*=( const CoefficientType& argument );
    
    /**
     *  @returns the number of non vanishing entries.
    **/
    size_t number_non_vanishing_entries() const;
    
    /**
     *  @returns true iff the vector is zero.
    **/
    bool is_zero() const;
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a vector.
     *  In order to do so call resize(0);
     */ 
    void resize ( const size_t dimension, const bool = true );
    
    /**
     *  @returns the number of entries.
     */
    size_t size() const;

    /**
     *  Fills every entry with CoefficientT(0).
     */
    void clear();
 
    /**
     *  @returns the homology class of the vector (with respect to differentials).
    **/
    ThisType homology_class(
            const MatrixField< CoefficientType >& base_changes_kernel,
            const MatrixField< CoefficientType >& base_changes_image
        ) const;
    
    // grant std::ostream access in order to print matrices to ostreams.
    template< class T >
    friend std::ostream& operator<< ( std::ostream& stream, const VectorField<T> & vector );

protected:
    VectorStorageType data; ///< This realizes the data.
    size_t dim;    ///< The number of entries.
    
    // In order to save vectors we have to grad boost::serialization::access access.
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int ) ///< Implements the serialization.
    {
        ar & dim & data;
    }
};

template< class CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const VectorField< CoefficientT > & vector );

////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Class for a vector with boolean coefficients
 * A VectorBool has the same structure as the VectorField,
 * but for a better performance, here each row of the vector is a dynamic_bitset.
 */
class VectorBool
{
public:
    typedef bool CoefficientType;
    typedef VectorBool ThisType;
    typedef boost::dynamic_bitset<> VectorStorageType;    ///< This realizes the implementation of the data.

    /**
     *  Creates a vector of dimension zero.
     */
    VectorBool();
    
    /**
     *  Creates a vector of a given dimension.
     *  In order to get a zero vector you may use clear();
     */
    VectorBool( const size_t dimension );

    /**
     *  In order to read elements of the vector \f$ v \f$ you want to use this function.
     *  @return The function returns a copy of the given entry.
     *  @todo throw an exception if necessary i.e. if \f$ v_i\f$ is not a valid entry.
     */
    bool operator()( const size_t i );
    
    /**
     *  In order to keep constness, we go the usual way and implement the function at.
     *  You may want to take a look at the at()-methods of the standard containers like std::vector.
     */
    bool at( const size_t i ) const;

    /**
     * Adds 1 to the entry i of the vector.
     * @note We need this function for writing access to the vector
     * since the dynamic bitsets does not provide a reasonable reference operator.
     */
    void add_entry( const size_t i );

    /**
     *  adds a given vector to this vector.
     */
    ThisType& operator+=( const ThisType& argument );
    
    /**
     *  substracts a given vector from this vector.
     */
    ThisType& operator-=( const ThisType& argument );
    
    /**
     *  multiply this vector with a given scalar.
     */
    ThisType& operator*=( const CoefficientType& argument );
    
    /**
     *  As our implementation mimes ublas::matrix we use the same (awkward) method to delete a vector.
     *  In order to do so call resize(0);
     */
    void resize ( const size_t dimension, const bool = true );

    /**
     *  @returns the number of entries.
     */
    size_t size() const;
    
    /**
     *  Fills every entry with 0.
     */
    void clear();
    
    // grant std::ostream access in order to print vectors to ostreams.
    friend std::ostream& operator<< ( std::ostream& stream, const VectorBool & vector);

protected:
    VectorStorageType data; ///< This realizes the data.
    size_t dim;             ///< The number of entries.

    // In order to save matrices we have to grad boost::serialization::access access.
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int ) ///< Implements the serialization.
    {
        ar & dim & data;
    }
};

std::ostream& operator<< ( std::ostream& stream, const VectorBool & vector);

////////////////////////////////////////////////////////////////////////////////////////////////////

template< class MatrixT, class VectorT >
void apply_base_changes_kernel( const MatrixT& m, VectorT& v );

template<>
void apply_base_changes_kernel( const MatrixBool& m, VectorBool& v );

template< class MatrixT, class VectorT >
void apply_base_changes_image( const MatrixT& m, VectorT& v );

template<>
void apply_base_changes_image( const MatrixBool& m, VectorBool& v );

template< class MatrixT, class VectorT >
std::vector< VectorT > compute_base_of_kernel( const MatrixT& m );

template<>
std::vector< VectorBool > compute_base_of_kernel( const MatrixBool& m );

template< class MatrixT, class VectorT >
VectorT matrix_vector_product( const MatrixT& m, const VectorT& v );

template< class MatrixT, class VectorT >
bool matrix_vector_product_vanishes( const MatrixT& m, const VectorT& v );

#endif // VECTOR_FIELD_HPP

