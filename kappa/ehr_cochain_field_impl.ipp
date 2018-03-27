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


#include "ehr_cochain_field.hpp"

// Delegate Constructur
template< typename CoefficientT >
EhrCochainField< CoefficientT >::EhrCochainField( const uint32_t genus, const uint32_t num_punct, const uint32_t cohom_deg ) :
        EhrCochainField< CoefficientT >(genus, num_punct, cohom_deg, false, "")
{
}

// Delegate Constructur
template< typename CoefficientT >
EhrCochainField< CoefficientT >::EhrCochainField( const uint32_t genus, const uint32_t num_punct, const uint32_t cohom_deg, const bool radial_model_used ) :
        EhrCochainField< CoefficientT >(genus, num_punct, cohom_deg, radial_model_used, "")
{
}

template< typename CoefficientT >
EhrCochainField< CoefficientT >::EhrCochainField( const uint32_t genus, const uint32_t num_punct, const uint32_t cohom_deg, const bool radial_model_used, const std::string& the_name ) :
        VectorType(), g(genus), m(num_punct), p(cohom_deg), radial(radial_model_used), name(the_name)
{
    basis = load_ehr_basis<SymGrpTuple>(g, m, p, radial);
    VectorType::resize( basis.size() );
}

template< typename CoefficientT >
typename EhrCochainField< CoefficientT >::CoefficientType & EhrCochainField< CoefficientT >::operator()( const SymGrpTuple& t )
{
    const auto res = basis.id_of(t);
    if( res == -1 )
    {
        std::cout << "Error: " << t << " is no basis element." << std::endl;
    }
    return VectorType::operator()( res );
}

template< typename CoefficientT >
const typename EhrCochainField< CoefficientT >::CoefficientType& EhrCochainField< CoefficientT > :: at( const SymGrpTuple& t ) const
{
    return VectorType::at( basis.id_of(t) );
}

template< typename CoefficientT >
std::string EhrCochainField< CoefficientT >::set_name( const std::string& new_name )
{
    return name = new_name;
}

template< typename CoefficientT >
void EhrCochainField< CoefficientT >::add_kappa_dual( const CoefficientType& c, const SymGrpTuple& t )
{
    VectorType::operator+=( kappa_dual< VectorType >( c, t, basis ) );
}

template< typename CoefficientT >
uint32_t EhrCochainField< CoefficientT >::get_g() const
{
    return g;
}

template< typename CoefficientT >
uint32_t EhrCochainField< CoefficientT >::get_m() const
{
    return m;
}

template< typename CoefficientT >
uint32_t EhrCochainField< CoefficientT >::get_p() const
{
    return p;
}

template< typename CoefficientT >
bool EhrCochainField< CoefficientT >::get_radial() const
{
    return radial;
}

template< typename CoefficientT >
std::string EhrCochainField< CoefficientT >::get_name() const
{
    return name;
}

template< typename CoefficientT >
const EhrBasis<SymGrpTuple>& EhrCochainField< CoefficientT >::get_basis_reference() const
{
    return basis;
}

template< typename CoefficientT >
EhrCochainField< CoefficientT > operator*( const EhrCochainField< CoefficientT >&  x, const EhrCochainField< CoefficientT >& y )
{
    typedef EhrCochainField< CoefficientT > CochainType;
    CochainType res( x.get_g() + y.get_g(), x.get_m() + y.get_m(), x.get_p() + y.get_p() );

    res.set_name( x.get_name() + y.get_name() );

    for( const auto cell_x : x.get_basis_reference().basis )
    {
        const CoefficientT& coeff_x = x.at( cell_x );
        if( coeff_x != CoefficientT(0) )
        {
            for( const auto cell_y : y.get_basis_reference().basis )
            {
                const CoefficientT& coeff_y = y.at( cell_y );
                if( coeff_y != CoefficientT(0) )
                {
                    res( cell_x * cell_y ) = coeff_x * coeff_y;
                }
            }
        }
    }

    return res;
}

template< typename CoefficientT >
std::ostream& operator<< ( std::ostream& stream, const EhrCochainField< CoefficientT > & cochain )
{
    return stream
        << "name = " << cochain.name
        << ", g = " << cochain.g
        << ", m = " << cochain.m
        << ", p = " << cochain.p
        << ", " << (cochain.radial == true ? "radial" : "parallel" ) << " version"
        << ", representing vector = "
        << static_cast< const typename EhrCochainField< CoefficientT >::VectorType & >(cochain);
}
