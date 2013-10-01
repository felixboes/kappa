#include "matrix_fp.hpp"

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 * Description:  initialize static members
 *--------------------------------------------------------------------------------------
 */
int Fp::base;
std::vector<int> Fp::inv;

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: set_modulus
 * Description:  calculate the inversetable of the ring Z/(p^k)
 *--------------------------------------------------------------------------------------
 */

void Fp::set_modulus(unsigned int p, unsigned int k)
{
    // Wir benutzen: Der Koeffizient i ist in Z/(p^k) genau dann invertierbar, wenn i%p != 0 ist.
    // Nach dem Satz von Euler ist dann a^{p^k - p^{k-1} - 1} = a^{-1} in Z/(p^k)
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

void Fp::clean_up()
{
    base = 0;
    inv.clear();
}

Fp::Fp(const int m) : n( ((m%base)+base)%base ) {}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: print_modulus()
 * Description:  print out the moduls
 *--------------------------------------------------------------------------------------
 */

void const Fp::print_modulus()
{
    std::cout << "The Coefficientring is F_" << base << ".\n";
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: printf_inversetable()
 * Description:  print out the inversetable
 *--------------------------------------------------------------------------------------
 */

void const Fp::print_inversetable()
{
    std::cout << "Printing inversetable:\n";
    for(int i = 1; i < base; i++)
    std::cout << i  << " " << inv[i] << "\n";
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: is_invertible
 * Description:  true iff n is invertible in Z mod p^k
 *--------------------------------------------------------------------------------------
 */

bool const Fp::is_invertible()
{
    return inv[((n % base)+base)%base];
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: inverse
 * Description:  returns the inverse of a Fpicient
 *--------------------------------------------------------------------------------------
 */

Fp const Fp::inverse()
{
    return Fp(((n % base)+base)%base);
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: print_Fp
 * Description:     prints our a Fpicient
 *--------------------------------------------------------------------------------------
 */
void const Fp::show()
{
    std::cout << ((n%base)+base)%base;
}

/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: print_Fp
 * Description:  prints a Fp to a String
 *--------------------------------------------------------------------------------------
 */
std::string Fp::get_string()
{
    std::ostringstream tmp;
    tmp << ((n%base)+base)%base;
    return tmp.str();
}


/*
 *--------------------------------------------------------------------------------------
 *       Class:  Fp
 *      Method:  Fp :: operator
 * Description:     operator+, ect...
 *--------------------------------------------------------------------------------------
 */

bool Fp::operator==(const int b) const
{
    return (n - b) % base == 0;
}

Fp& Fp::operator=(const int a)
{
    n = ((a%base)+base)%base;
    return *this;
}

Fp& Fp::operator +=(Fp a)
{
    n = (((n+a.n)%base)+base)%base;
    return *this;
}

Fp& Fp::operator -=(Fp a)
{
    n = (((n-a.n)%base)+base)%base;
    return *this;
}

Fp& Fp::operator *=(Fp a)
{
    n = (((n*a.n)%base)+base)%base;
    return *this;
}

Fp Fp::operator-() const
{
    Fp r(-n);
    return r;
}

Fp::operator bool() const
{
    return n%base;
}


Fp operator+(Fp a, Fp b)
{
    Fp r(a);
    return r += b;
}

Fp operator-(Fp a, Fp b)
{
    Fp r(a);
    return r -= b;
}

Fp operator*(Fp a, Fp b)
{
    Fp r(a);
    return r *= b;
}

Fp operator*(Fp a, const int b)
{
    Fp r(a);
    return r *= Fp(b);
}
