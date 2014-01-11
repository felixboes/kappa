#ifndef FIELD_COEFFICIENTS_HPP
#define FIELD_COEFFICIENTS_HPP

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
    Zm(const int8_t m = 0);
    static void set_modulus(const uint8_t prime, const uint8_t expo = 1);
    static void const print_modulus();
    static void const print_inversetable();
    bool const is_invertible();
    Zm const inverse();
    void const show();
    static void clean_up();
    static bool is_field();
    
    // arithmetic operators
    bool operator==(const int8_t) const;
    bool operator==(const Zm a) const;
    Zm& operator=(const int8_t);
    Zm& operator+=(const Zm);
    Zm& operator-=(const Zm);
    Zm& operator*=(const Zm);
    Zm& operator/=(const Zm);
    Zm operator-() const;
    operator bool() const;  // false iff the coefficient is zero
    
    friend std::ostream& operator<< (std::ostream& stream, const Zm& coeff);
    
private:
    int8_t n;
    static uint8_t prim;
    static uint8_t expo;
    static int8_t base;        // base =  p^k
    static std::vector<int8_t> inv;
    operator int() const;
    operator unsigned() const;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization of a coefficient.
    {
        ar & n;
    }
};

bool operator !=( const Zm, const Zm );
Zm operator+(const Zm, const Zm);
Zm operator-(const Zm, const Zm);
Zm operator*(const Zm, const Zm);
Zm operator/(const Zm, const Zm);

Zm operator*(const Zm, const int8_t);



#endif // ENDIF FIELD_COEFFICIENTS_HPP
