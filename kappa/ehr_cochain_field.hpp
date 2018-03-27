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


#ifndef EHR_COCHAINFIELD_H
#define EHR_COCHAINFIELD_H

#include <libhomology/homology.hpp>

#include "sym_grp_tuple.hpp"
#include "ehr_basis.hpp"
#include "operationtester.hpp"

template< class CoefficientT >
class EhrCochainField : public VectorField< CoefficientT >
{
public:
    typedef CoefficientT CoefficientType;
    typedef EhrCochainField< CoefficientType > ThisType;
    typedef typename VectorField< CoefficientType > :: ThisType VectorType;
    typedef typename VectorField< CoefficientType > :: VectorStorageType VectorStorageType;

    EhrCochainField( const uint32_t genus, const uint32_t num_punct, const uint32_t cohom_deg );
    EhrCochainField( const uint32_t genus, const uint32_t num_punct, const uint32_t cohom_deg, const bool radial_model_used );
    EhrCochainField( const uint32_t genus, const uint32_t num_punct, const uint32_t cohom_deg, const bool radial_model_used, const std::string& the_name );

    CoefficientType & operator()( const SymGrpTuple& t );

    const CoefficientType& at( const SymGrpTuple& t ) const;

    std::string set_name( const std::string& new_name );

    void add_kappa_dual( const CoefficientType& c, const SymGrpTuple& t );

    template< class T >
    friend ThisType operator*( const ThisType&, const ThisType& );

    uint32_t get_g() const;

    uint32_t get_m() const;

    uint32_t get_p() const;

    bool get_radial() const;

    std::string get_name() const;

    const EhrBasis& get_basis_reference() const;

    // grant std::ostream access in order to print matrices to ostreams.
    template< class T >
    friend std::ostream& operator<< ( std::ostream& stream, const EhrCochainField<T> & cochain );

protected:
    EhrBasis basis;
    uint32_t g;
    uint32_t m;
    uint32_t p;
    bool radial;
    std::string name;
};

template< typename CoefficientT >
EhrCochainField< CoefficientT > operator*( const EhrCochainField< CoefficientT >&  x, const EhrCochainField< CoefficientT >& y );

#endif