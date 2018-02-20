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


#include "field_coefficients.hpp"
#include "field_coefficients_impl.ipp"

#include "homology.hpp"

/* Force template instantiation for used types */

template class ZmBase<>;
template std::ostream& operator<< (std::ostream& stream, const ZmBase<>& coeff);

template bool operator !=( const ZmBase<>&, const ZmBase<>& );
template bool operator !=( const ZmBase<>&, const ZmBase<>::BaseType );
template ZmBase<> operator+(const ZmBase<>&, const ZmBase<>&);
template ZmBase<> operator-(const ZmBase<>&, const ZmBase<>&);
template ZmBase<> operator*(const ZmBase<>&, const ZmBase<>&);
template ZmBase<> operator/(const ZmBase<>&, const ZmBase<>&);
template ZmBase<> operator*(const ZmBase<>&, const ZmBase<>::BaseType);

template ZmBase<> di (const ZmBase<>&, const ZmBase<>&);
template ZmBase<> mod(const ZmBase<>&, const ZmBase<>&);
template ZmBase<> gcd(const ZmBase<>&, const ZmBase<>&);
template std::pair<ZmBase<>, ZmBase<> >bezout(const ZmBase<>&, const ZmBase<>&);
