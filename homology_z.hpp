#ifndef HOMOLOGY_Z_HPP
#define HOMOLOGY_Z_HPP

#include <stdint.h>
#include <map>

/**
 TODO: Check the torsion part!

 *  If we work over \f$ \mathbb{Z} \f$, a homology module \f$ H \f$ is a finite sum of free and torsion modules.
 *  \f[
 *      H = \mathbb{Z}^N \oplus \mathbb{Z}/a_1\mathbb{Z} \oplus \ldots \oplus \mathbb{Z}/a_r\mathbb{Z}
 *  \f]
 *  At the moment we remember the numbers \f$N\f$ and the \f$a_i\f$.
 *  A chosen basis or base change is not stored.
 */
template< class Coefficient > class HomologyModuleZ
{
public:

    /// trivial constructor
    HomologyModuleZ();

    /// Access the dimension of the free part.
    uint32_t free_part();

    /// Access the number of torsion modules.
    uint32_t tors_part();

    /// Access the \f$n\f$-th torsion coefficient \f$ a_n \f$.
    Coefficient &operator[] (uint32_t n);

private:

    /// The number \f$N\f$
    uint32_t free_part;

    /// The number \f$r\f$
    uint32_t tors_part;

    std::map< uint32_t, Coefficient > tors_coeff;
};

template< class Coefficient > class HomologyComplexZ
{
public:
    HomologyComplexZ();
    HomologyModuleZ &operator[] ( uint32_t n ) { return module[n]; }

private:
    std::map< uint32_t, HomologyModuleZ > module;
};

#endif // HOMOLOGY_Z_HPP
