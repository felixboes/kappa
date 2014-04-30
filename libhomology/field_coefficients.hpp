#ifndef FIELD_COEFFICIENTS_HPP
#define FIELD_COEFFICIENTS_HPP

// Description:
//
// This header defines the coefficient rings Q and Z/mZ.
// Befor using Z/mZ coefficients you must use the static directive set_modulus.

#include <cinttypes>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <vector>

#include <boost/serialization/access.hpp>

#include "parallelization.hpp"

/**
 *  The coefficient ring \f$ \mathbb{Q} \f$.
 */
typedef mpq_class Q;

/**
 *  The coefficient ring \f$ \mathbb{Z} / m \mathbb{Z} \f$.
 */
class Zm{
public:
    Zm(const int8_t m = 0);    ///< The default constructor creates a coefficient with value zero.
    static void set_modulus(const uint8_t prime, const uint8_t expo = 1);    ///< Befor using Zm coefficients, you must define the modulus i.e. m = p^e.
    static void const print_modulus();    ///< Print the modulus to std::cout.
    static void const print_inversetable();    ///< Print the table of invertible elements to std::cout.
    bool const is_invertible(); ///< @returns true iff the coefficient is invertible.
    Zm const inverse();         ///< @returns the inverse of a given coefficient. If the coefficient is not invertible the value zero is returned.
    static void clean_up();     ///< Clean up all static data e.g. the table of invertible elements.
    static bool is_field();     ///< @returns true iff m = p^e is a prime number. @todo: primeness of p is not yet verified.
    
    // arithmetic operators
    bool operator==(const int8_t) const;    ///< compare a Zm with an int8_t
    bool operator==(const Zm a) const;      ///< compare a Zm with another Zm.
    Zm& operator=(const int8_t);    ///< Assignement.
    Zm& operator+=(const Zm);       ///< Adding.
    Zm& operator-=(const Zm);       ///< Subtracting.
    Zm& operator*=(const Zm);       ///< Multiplying.
    Zm& operator/=(const Zm);       ///< Dividing. @todo: throw exeption if necessary.
    Zm operator-() const;           ///< Inverting additively.
    operator bool() const;          ///< @returns false iff the coefficient is zero.
    
    /**
     *  grant std::ostream access in order to print coefficients to ostreams like 'std::cout << Zm(44) << std::endl;'
     */
    friend std::ostream& operator<< (std::ostream& stream, const Zm& coeff);
    
private:
    int8_t n;   ///< This integer stores a representative of the residue class \f$ c = [n] in \mathbb{Z}/m\mathbb{Z} \f$.
    static uint8_t prim;    ///< We store the number m = p^e for all coeffients at once. Therefore we have to use set_modulus befor working with such coefficients.
    static uint8_t expo;    ///< We store the number m = p^e for all coeffients at once. Therefore we have to use set_modulus befor working with such coefficients.
    static int8_t base;     ///< m = base = p^k.
    static std::vector<int8_t> inv; ///< This vector stores the table of inverse elements.
    operator int() const;       ///< In order to cast a Zm coefficient \f$c\f$ to an integer we pick a representative \f$ 0 \le c < base\f$.
    operator unsigned() const;  ///< In order to cast a Zm coefficient \f$c\f$ to an unsigned integer we pick a representative \f$ 0 \le c < base\f$.
    
    /**
     *  In order to save Zm coefficients we have to grad boost::serialization::access access.
     */
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization of a coefficient.
    {
        ar & n;
    }
};

bool operator !=( const Zm, const Zm ); ///< Compare two Zm coeffiencts and return true iff they are different.
Zm operator+(const Zm, const Zm);   ///< Add two Zm Coeffients and return the result.
Zm operator-(const Zm, const Zm);   ///< Substract two Zm Coeffients and return the result.
Zm operator*(const Zm, const Zm);   ///< Multiply two Zm Coeffients and return the result.
Zm operator/(const Zm, const Zm);   ///< Divide two Zm Coeffients and return the result. @todo throw exception if necessary.

Zm operator*(const Zm, const int8_t);   ///< Multiply a Zm coefficient and an integer and return the result.

#endif // ENDIF FIELD_COEFFICIENTS_HPP
