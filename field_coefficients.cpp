#include "field_coefficients.hpp"

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 * Description:  initialize static members
 *--------------------------------------------------------------------------------------
 */
uint8_t Zm::prim;
uint8_t Zm::expo;
int8_t Zm::base;
std::vector<int8_t> Zm::inv;

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: set_modulus
 * Description:  calculate the inversetable of the ring Z/(p^k)
 *--------------------------------------------------------------------------------------
 */

void Zm::set_modulus(const uint8_t p, const uint8_t k)
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
    for(int8_t i = 2; i < base; i++)    // calculate modulo exponent
    {
        if(i % p != 0)
        {
            int32_t ex = m;
            int8_t factor = i;
            int8_t result = 1;
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

void Zm::clean_up()
{
    base = 0;
    inv.clear();
}

bool Zm::is_field()
{
    // todo check primeness.
    return (prim > 0 && expo == 1);
}

Zm::Zm(const int8_t m) : n( ((m%base)+base)%base ) {}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: print_modulus()
 * Description:  print out the moduls
 *--------------------------------------------------------------------------------------
 */

void const Zm::print_modulus()
{
    std::cout << "The Coefficientring is F_" << (int32_t)base << ".\n";
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: printf_inversetable()
 * Description:  print out the inversetable
 *--------------------------------------------------------------------------------------
 */

void const Zm::print_inversetable()
{
    std::cout << "Printing inversetable:\n";
    for(size_t i = 1; i < base; i++)
    std::cout << (int32_t)i  << " " << inv[i] << "\n";
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: is_invertible
 * Description:  true iff n is invertible in Z mod p^k
 *--------------------------------------------------------------------------------------
 */

bool const Zm::is_invertible()
{
    return inv[((n % base)+base)%base];
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: inverse
 * Description:  returns the inverse of a Coefficient
 *--------------------------------------------------------------------------------------
 */

Zm const Zm::inverse()
{
    return Zm( inv[(((n % base)+base)%base)]);
}

std::ostream& operator<< (std::ostream& stream, const Zm& coeff)
{
    return stream << (int32_t)((coeff.n %coeff.base)+coeff.base)%coeff.base;
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: operator
 * Description:     operator+, ect...
 *--------------------------------------------------------------------------------------
 */

bool Zm::operator==(const int8_t b) const
{
    return (n - b) % base == 0;
}

bool Zm::operator== (const Zm a) const
{
    return (n - a.n) % base == 0;
}

bool operator !=( const Zm a, const Zm b )
{
    return !( a.operator==(b) );
}


Zm& Zm::operator=(const int8_t a)
{
    n = ((a%base)+base)%base;
    return *this;
}

Zm& Zm::operator +=(Zm a)
{
    n = (((n+a.n)%base)+base)%base;
    return *this;
}

Zm& Zm::operator -=(Zm a)
{
    n = (((n-a.n)%base)+base)%base;
    return *this;
}

Zm& Zm::operator *=(Zm a)
{
    n = (((n*a.n)%base)+base)%base;
    return *this;
}

Zm& Zm::operator /=(Zm a)
{
    this->operator*= ( a.inverse() );
    return *this;
}

Zm Zm::operator-() const
{
    Zm r(-n);
    return r;
}

Zm::operator bool() const
{
    return n%base;
}

Zm operator+(Zm a, Zm b)
{
    Zm r(a);
    return r += b;
}

Zm operator-(Zm a, Zm b)
{
    Zm r(a);
    return r -= b;
}

Zm operator*(Zm a, Zm b)
{
    Zm r(a);
    return r *= b;
}

Zm operator/ (Zm a, Zm b)
{
    Zm r(a);
    return r /= b;
}

Zm operator*(Zm a, const int8_t b)
{
    Zm r(a);
    return r *= Zm(b);
}
