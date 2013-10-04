#ifndef MATRIX_ZM_HPP
#define MATRIX_ZM_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

/**
 *  The coefficient ring \f$ \mathbb{Z} / m \mathbb{Z} \f$.
 */
class Zm{
public:
    Zm(const int m = 0);
    static void set_modulus(const uint prime, const uint expo);
    static void const print_modulus();
    static void const print_inversetable();
    bool const is_invertible();
    Zm const inverse();
    void const show();
    std::string get_string();
    static void clean_up();
    static bool is_field();
    
    // arithmetic operators
    bool operator==(const int) const;
    Zm& operator=(const int);
    Zm& operator+=(const Zm);
    Zm& operator-=(const Zm);
    Zm& operator*=(const Zm);
    Zm operator-() const;
    operator bool() const;  // false iff the coefficient is zero
    
    friend std::ostream& operator<< (std::ostream& stream, const Zm& coeff);
    
private:
    int n;
    static uint prim;
    static uint expo;
    static int base;        // base =  p^k
    static std::vector<int> inv;
    operator int() const;
    operator unsigned() const;
};

Zm operator+(const Zm, const Zm);
Zm operator-(const Zm, const Zm);
Zm operator*(const Zm, const Zm);
Zm operator*(const Zm, const int);

std::ostream& operator<< (std::ostream& stream, const Zm& coeff);

typedef boost::numeric::ublas::matrix< Zm > MatrixZm;
typedef boost::numeric::ublas::matrix_slice< MatrixZm > MatrixZmSlice;
typedef boost::numeric::ublas::slice slice_zm;

#endif // MATRIX_Zm_HPP
