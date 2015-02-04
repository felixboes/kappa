#ifndef FIELD_COEFFICIENTS_HPP
#define FIELD_COEFFICIENTS_HPP

// Description:
//
// This header defines the coefficient rings Q and Z/mZ.
// Befor using Z/mZ coefficients you must use the static directive set_modulus.

#include <cstddef> // has to be included before gmpxx with gcc 4.9
#include <cinttypes>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#include <boost/serialization/access.hpp>

/**
 *  The coefficient ring \f$ \mathbb Q \f$.
 */
typedef mpq_class Q;

/**
 *  The coefficient ring \f$ \mathbb Z / m \mathbb Z \f$.
 */
template < typename base_type = int8_t >
class ZmBase{
public:
    typedef base_type BaseType;
    typedef ZmBase< BaseType > ThisType;
    
    ZmBase(const BaseType m = 0);    ///< The default constructor creates a coefficient with value zero.
    static void set_modulus(const BaseType prime, const BaseType expo = 1);    ///< Befor using ZmBase coefficients, you must define the modulus i.e. m = p^e.
    static void print_modulus();    ///< Print the modulus to std::cout.
    static void print_inversetable();    ///< Print the table of invertible elements to std::cout.
    static BaseType get_modulus();   ///< @returns the modulus.
    bool is_invertible() const; ///< @returns true iff the coefficient is invertible.
    ThisType inverse() const;         ///< @returns the inverse of a given coefficient. If the coefficient is not invertible the value zero is returned.
    static void clean_up();     ///< Clean up all static data e.g. the table of invertible elements.
    static bool is_field();     ///< @returns true iff m = p^e is a prime number. @todo: primeness of p is not yet verified.
    
    // arithmetic operators
    bool operator==(const BaseType) const;  ///< compare a ZmBase with an int8_t
    bool operator==(const ThisType) const;  ///< compare a ZmBase with another ZmBase.
    ThisType& operator= (const BaseType);   ///< Assignement.
    ThisType& operator+=(const ThisType);   ///< Adding.
    ThisType& operator-=(const ThisType);   ///< Subtracting.
    ThisType& operator*=(const ThisType);   ///< Multiplying.
    ThisType& operator/=(const ThisType);   ///< Dividing. @todo: throw exeption if necessary.
    ThisType operator-() const;             ///< Inverting additively.
    operator bool() const;          ///< @returns false iff the coefficient is zero.
    
    ThisType& di (const ThisType);
    ThisType& mod(const ThisType);
    
    /**
     *  grant std::ostream access in order to print coefficients to ostreams like 'std::cout << ZmBase(44) << std::endl;'
     */
    template < typename T >
    friend std::ostream& operator<< (std::ostream& stream, const ZmBase<T>& coeff);
    
protected:
    BaseType n;    ///< This integer stores a representative of the residue class \f$ c = [n] in \mathbb{Z}/m\mathbb{Z} \f$.
    static BaseType prim;    ///< We store the number m = p^e for all coeffients at once. Therefore we have to use set_modulus befor working with such coefficients.
    static BaseType expo;    ///< We store the number m = p^e for all coeffients at once. Therefore we have to use set_modulus befor working with such coefficients.
    static BaseType base;  ///< m = base = p^k.
    static std::vector<BaseType> inv; ///< This vector stores the table of inverse elements.
    operator int() const;       ///< In order to cast a ZmBase coefficient \f$c\f$ to an integer we pick a representative \f$ 0 \le c < base\f$.
    operator unsigned() const;  ///< In order to cast a ZmBase coefficient \f$c\f$ to an unsigned integer we pick a representative \f$ 0 \le c < base\f$.
    
    // In order to save ZmBase coefficients we have to grad boost::serialization::access access.
    friend class boost::serialization::access;
    
    template < class Archive >
    void serialize( Archive &ar, const unsigned int ) ///< Implements the serialization of a coefficient.
    {
        ar & n;
    }
};

template < typename base_type = int8_t > bool operator !=( const ZmBase<base_type>, const ZmBase<base_type> ); ///< Compare two ZmBase coeffiencts and return true iff they are different.
template < typename base_type = int8_t > ZmBase<base_type> operator+(const ZmBase<base_type>, const ZmBase<base_type>);   ///< Add two ZmBase Coeffients and return the result.
template < typename base_type = int8_t > ZmBase<base_type> operator-(const ZmBase<base_type>, const ZmBase<base_type>);   ///< Substract two ZmBase Coeffients and return the result.
template < typename base_type = int8_t > ZmBase<base_type> operator*(const ZmBase<base_type>, const ZmBase<base_type>);   ///< Multiply two ZmBase Coeffients and return the result.
template < typename base_type = int8_t > ZmBase<base_type> operator/(const ZmBase<base_type>, const ZmBase<base_type>);   ///< Divide two ZmBase Coeffients and return the result. @todo throw exception if necessary.
template < typename base_type = int8_t > ZmBase<base_type> operator*(const ZmBase<base_type>, const base_type);   ///< Multiply a ZmBase coefficient and an integer and return the result.

template < typename base_type = int8_t > ZmBase<base_type> di (const ZmBase<base_type>, const ZmBase<base_type>);
template < typename base_type = int8_t > ZmBase<base_type> mod(const ZmBase<base_type>, const ZmBase<base_type>);
template < typename base_type = int8_t > ZmBase<base_type> gcd(const ZmBase<base_type>, const ZmBase<base_type>);
template < typename base_type = int8_t > std::pair<ZmBase<base_type>, ZmBase<base_type> >bezout(const ZmBase<base_type>, const ZmBase<base_type>);

typedef ZmBase<> Zm;

#include "field_coefficients.ipp"

#endif // ENDIF FIELD_COEFFICIENTS_HPP
