#include "matrix_zm.hpp"

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 * Description:  initialize static members
 *--------------------------------------------------------------------------------------
 */
uint32_t Zm::prim;
uint32_t Zm::expo;
int32_t Zm::base;
std::vector<int> Zm::inv;

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: set_modulus
 * Description:  calculate the inversetable of the ring Z/(p^k)
 *--------------------------------------------------------------------------------------
 */

void Zm::set_modulus(const uint p, const uint k)
{
    // Wir benutzen: Der Koeffizient i ist in Z/(p^k) genau dann invertierbar, wenn i%p != 0 ist.
    // Nach dem Satz von Euler ist dann a^{p^k - p^{k-1} - 1} = a^{-1} in Z/(p^k)
    prim = p;
    expo = k;
    base = pow(p,k);
    inv.clear();
    inv.push_back(0);
    inv.push_back(1);
    int m = pow(p,k) - pow(p,k-1) - 1;
    for(int i = 2; i < base; i++)    // calculate modulo exponent
    {
        if(i % p != 0)
        {
            int ex = m;
            int factor = i;
            int result = 1;
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
    return (prim > 0 && expo == 1);
}

Zm::Zm(const int m) : n( ((m%base)+base)%base ) {}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: print_modulus()
 * Description:  print out the moduls
 *--------------------------------------------------------------------------------------
 */

void const Zm::print_modulus()
{
    std::cout << "The Coefficientring is F_" << base << ".\n";
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
    for(int i = 1; i < base; i++)
    std::cout << i  << " " << inv[i] << "\n";
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
 * Description:  returns the inverse of a Zmicient
 *--------------------------------------------------------------------------------------
 */

Zm const Zm::inverse()
{
    return Zm(((n % base)+base)%base);
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: print_Zm
 * Description:     prints our a Zmicient
 *--------------------------------------------------------------------------------------
 */
void const Zm::show()
{
    std::cout << ((n%base)+base)%base;
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: print_Zm
 * Description:  prints a Zm to a String
 *--------------------------------------------------------------------------------------
 */
std::string Zm::get_string()
{
    std::ostringstream tmp;
    tmp << ((n%base)+base)%base;
    return tmp.str();
}

std::ostream& operator<< (std::ostream& stream, const Zm& coeff)
{
    return stream << coeff.n;
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Zm
 *      Method:  Zm :: operator
 * Description:     operator+, ect...
 *--------------------------------------------------------------------------------------
 */

bool Zm::operator==(const int b) const
{
    return (n - b) % base == 0;
}

Zm& Zm::operator=(const int a)
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

Zm operator*(Zm a, const int b)
{
    Zm r(a);
    return r *= Zm(b);
}
