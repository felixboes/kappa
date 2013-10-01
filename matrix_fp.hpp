#ifndef MATRIX_FP_HPP
#define MATRIX_FP_HPP

#include <boost/numeric/ublas/matrix.hpp>

class Fp{
public:
    Fp(const int m = 0) : n(m) {}
    static void set_modulus(const unsigned prime, const unsigned int expo);
    static void const print_modulus();
    static void const print_inversetable();
    bool const is_invertible();
    Fp const inverse();
    void const show();
    std::string getString();
    static void clean_up();

    // arithmetic operators
    bool operator==(const int) const;
    Fp& operator=(const int);
    Fp& operator+=(const Fp);
    Fp& operator-=(const Fp);
    Fp& operator*=(const Fp);
    Fp operator-() const;
    operator bool() const;  // false iff the coefficient is zero

private:
    int n;
    static int base;        // base =  p^k
    static std::vector<int> inv;
    operator int() const;
    operator unsigned() const;
};

Fp operator+(const Fp, const Fp);
Fp operator-(const Fp, const Fp);
Fp operator*(const Fp, const Fp);
Fp operator*(const Fp, const int);

typedef boost::numeric::ublas::matrix< Fp > MatrixFp;

#endif // MATRIX_FP_HPP
