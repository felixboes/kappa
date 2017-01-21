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
