#include "field_coefficients.hpp"

// Initialize static memebers
template < typename base_type >
uint8_t ZmBase<base_type>::prim;
template < typename base_type >
uint8_t ZmBase<base_type>::expo;
template < typename base_type >
typename ZmBase<base_type>::BaseType ZmBase<base_type>::base;
template < typename base_type >
typename std::vector< typename ZmBase<base_type>::BaseType > ZmBase<base_type>::inv;

// Implementation of (member-) functions
template < typename base_type >
void ZmBase<base_type>::set_modulus(const uint8_t p, const uint8_t k)
{
    // A coefficient Z/(p^k) is invertible iff i%p != 0 .
    // By Eulers theorem a^{p^k - p^{k-1} - 1} = a^{-1} in Z/(p^k)
    prim = p;
    expo = k;
    base = pow(p,k);
    inv.clear();
    inv.push_back(0);   // Here the inverse of zero is zero.
    inv.push_back(1);   // The inverse of 1 is cleary 1.
    
    int32_t m = pow(p,k) - pow(p,k-1) - 1;
    for(BaseType i = 2; i < base; i++)    // calculate modulo exponent
    {
        if(i % p != 0)
        {
            int32_t ex = m;
            BaseType factor = i;
            BaseType result = 1;
            while(ex > 0)
            {
                if(ex % 2 == 1)
                {
                    result = (result * factor) % base;
                }

                factor = (factor * factor) % base;
                ex /= 2;
            }
            inv.push_back(result);
        }
        else
        {
            inv.push_back(0);
        }
    }
}

template < typename base_type >
typename ZmBase<base_type>::BaseType ZmBase<base_type>::get_modulus()
{
    return prim;
}

template < typename base_type >
void ZmBase<base_type>::clean_up()
{
    base = 0;
    inv.clear();
}

template < typename base_type >
bool ZmBase<base_type>::is_field()
{
    // todo check primeness.
    return (prim > 0 && expo == 1);
}

template < typename base_type >
ZmBase<base_type>::ZmBase(const BaseType m) : n( ((m%base)+base)%base ) {}

template < typename base_type >
void ZmBase<base_type>::print_modulus()
{
    std::cout << "The Coefficientring is F_" << (int64_t)base << std::endl;
}

template < typename base_type >
void ZmBase<base_type>::print_inversetable()
{
    std::cout << "Printing inversetable:" << std::endl;
    for(BaseType i = 1; i < base; i++)
    {
        std::cout << (int64_t)i  << " " << (int64_t)inv[i] << "\n";
    }
}

template < typename base_type >
bool ZmBase<base_type>::is_invertible() const
{
    return inv[((n % base)+base)%base];
}

template < typename base_type >
ZmBase<base_type> ZmBase<base_type>::inverse() const
{
    return ZmBase<BaseType>( inv[(((n % base)+base)%base)]);
}

template < typename base_type = uint8_t >
std::ostream& operator<< (std::ostream& stream, const ZmBase<base_type>& coeff)
{
    return stream << (int64_t)((coeff.n %coeff.base)+coeff.base)%coeff.base;
}

template < typename base_type >
bool ZmBase<base_type>::operator==(const BaseType b) const
{
    return (n - b) % base == 0;
}

template < typename base_type >
bool ZmBase<base_type>::operator== (const ThisType a) const
{
    return (n - a.n) % base == 0;
}

template < typename base_type >
bool operator !=( const ZmBase<base_type> a, const ZmBase<base_type> b )
{
    return !( a.operator==(b) );
}

template < typename base_type >
ZmBase<base_type>& ZmBase<base_type>::operator=(const BaseType a)
{
    n = ((a%base)+base)%base;
    return *this;
}

template < typename base_type >
ZmBase<base_type>& ZmBase<base_type>::operator +=(ThisType a)
{
    n = (((n+a.n)%base)+base)%base;
    return *this;
}

template < typename base_type >
ZmBase<base_type>& ZmBase<base_type>::operator -=(ThisType a)
{
    n = (((n-a.n)%base)+base)%base;
    return *this;
}

template < typename base_type >
ZmBase<base_type>& ZmBase<base_type>::operator *=(ThisType a)
{
    n = (((n*a.n)%base)+base)%base;
    return *this;
}

template < typename base_type >
ZmBase<base_type>& ZmBase<base_type>::operator /=(ThisType a)
{
    this->operator*= ( a.inverse() );
    return *this;
}

template < typename base_type >
ZmBase<base_type> ZmBase<base_type>::operator-() const
{
    ZmBase<BaseType> r(-n);
    return r;
}

template < typename base_type >
ZmBase<base_type>::operator bool() const
{
    return n%base;
}

template < typename base_type >
ZmBase<base_type> operator+(ZmBase<base_type> a, ZmBase<base_type> b)
{
    ZmBase<base_type> r(a);
    return r += b;
}

template < typename base_type >
ZmBase<base_type> operator-(ZmBase<base_type> a, ZmBase<base_type> b)
{
    ZmBase<base_type> r(a);
    return r -= b;
}

template < typename base_type >
ZmBase<base_type> operator*(ZmBase<base_type> a, ZmBase<base_type> b)
{
    ZmBase<base_type> r(a);
    return r *= b;
}

template < typename base_type >
ZmBase<base_type> operator/ (ZmBase<base_type> a, ZmBase<base_type> b)
{
    ZmBase<base_type> r(a);
    return r /= b;
}

template < typename base_type >
ZmBase<base_type> operator*(ZmBase<base_type> a, const base_type b)
{
    ZmBase<base_type> r(a);
    return r *= ZmBase<base_type>(b);
}