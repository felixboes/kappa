#ifndef MONOCOMPLEX_H
#define MONOCOMPLEX_H

#include <stdint.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <vector>

#include <homology>

#include "factorial.hpp"
#include "tuple.hpp"


/**
    The MonoBasis keeps track of the basis elements.
**/
struct MonoBasis
{
    /// Add a basis element.
    void add_basis_element (Tuple& t);
    
    /// output stream
    friend std::ostream& operator<< (std::ostream& stream, const MonoBasis& basis);
    
    /// Returns the number of basis elements.
    uint64_t size();

    /// Stores the orderd basis.
    std::vector< Tuple > basis;
};

std::ostream& operator<< (std::ostream& stream, const MonoBasis& basis);

/**
    Kettenkomplex der von den monotonen Tupeln von Transpositionen erzeugt wird.
**/
template< class MatrixComplex >
class MonoComplex
{
public:
    MonoComplex(uint32_t genus, uint32_t num_punctures);

    // void show();        ///< Gibt den Kettenkomplex auf der Standardausgabe aus.
    // void show_size();   ///< Gibt die Groesse der einzelnen Kettenmoduln aus.
    /**
        Siehe Erklaerung 1.
    **/
    void gen_bases(uint32_t l, uint32_t p, Tuple& tuple);
    void gen_differential(int32_t p);
    void gen_differentials();
    
    //void del2_kappa(uint32_t p); ///< Die Funktion \f$ kappa_p \f$, wobei p der horizontale Grad ist.

private:
    uint32_t g;     ///< Geschlecht
    uint32_t m;     ///< Anzahl der Punktierungen
    uint32_t h;     ///< h = 2*g+m
    
    MatrixComplex matrix_complex;
    std::map< uint32_t, MonoBasis > basis_complex;   ///< modul[n] ist der n-te Kettenmodul.
};

#endif // MONOCOMPLEX_H
