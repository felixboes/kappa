#ifndef MONOKOMPLEX_H
#define MONOKOMPLEX_H

#include <stdint.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <vector>

#include "factorial.hpp"
#include "tuple.hpp"

/**
    Ein frei erzeugter Modul.
    Die Basis sind die monotonen Tupeln von Transpositionen.
    Ein Tupel ist genau dann monoton, wenn die folgende Bedingung erfuellt ist:
    \f[
        ht( \tau_q ) \ge ht( \tau_{q-1} ) \ge \ldots \ge ht( \tau_1 ) \,.
    \f]
    Gespeichert wird nur die Basis.
**/
struct MonoModul
{
    void ergaenze (Tupel& erzeuger);     ///< Ergaenzt die Basis des Moduls um einen Erzeuger.
    void show();        ///< Gibt die Basis auf der Standardausgabe aus.
    uint64_t size();    ///< Anzahl der Erzeuger

    std::vector< Tupel > omega;     ///< Die Menge der monotonen Tupel
};


/**
    Kettenkomplex der von den monotonen Tupeln von Transpositionen erzeugt wird.
**/
class MonoKomplex
{
public:
    MonoKomplex(uint32_t _g, uint32_t _m, bool _externe);

    void show();        ///< Gibt den Kettenkomplex auf der Standardausgabe aus.
    void show_size();   ///< Gibt die Groesse der einzelnen Kettenmoduln aus.
    /**
        Siehe Erklaerung 1.
    **/
    void gen_modulen_extern(uint32_t l, uint32_t p, Tupel& anfangsstueck); // externe Version
    void gen_modulen_intern(uint32_t l, uint32_t p, Tupel& anfangsstueck, uint32_t intpkt ); // interne Version
    void del2_kappa(uint32_t p); ///< Die Funktion \f$ kappa_p \f$, wobei p der horizontale Grad ist.

private:
    uint32_t g;     ///< Geschlecht
    uint32_t m;     ///< Anzahl der Punktierungen
    uint32_t h;     ///< h = 2*g+m
    bool externe;   ///< Genau dann wahr, wenn wir externe Punktierungen benutzen.
    std::map< uint32_t, MonoModul > modul;   ///< modul[n] ist der n-te Kettenmodul.
};

#endif // MONOKOMPLEX_H
