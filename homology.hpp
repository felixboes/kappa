#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

#include <stdint.h>
#include <map>

template< class Coefficient > class HomologyModule
{
public:
    HomologyModule();

    Coefficient &operator[] (uint32_t n) { return tors_coeff[n]; }

    uint32_t free_part;
    uint32_t tors_part;

private:
    std::map< uint32_t, Coefficient > tors_coeff;
};

template< class Coefficient > class HomologyComplex
{
public:
    HomologyComplex();
    HomologyModule &operator[] ( uint32_t n ) { return module[n]; }

private:
    std::map< uint32_t, HomologyModule > module;
};

#endif // HOMOLOGY_HPP
