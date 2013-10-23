#ifndef TUPLE_HPP
#define TUPLE_HPP

#include <cinttypes>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

typedef std::pair< uint8_t, uint8_t > Transposition;

struct PermutationType
{
    PermutationType() : num_cycles(0), num_punctures(0) {}
    
    uint32_t num_cycles;
    uint32_t num_punctures;
};

/**
 *  A Tuple of Transpositions \f$ ( \tau_h \mid \ldots \mid \tau_1 ) \f$.
 *  
 *  Convention: We want an easy-to-check condition, to test if a Tuple is in a bad state.
 *  We require, that good Tuples use only the symbols \f$ 1 \f$ to \f$ p \f$.
**/
class Tuple
{
public:
    /**
     *  Construct a Tuple of norm h, which has to be filled.
     */
    Tuple(size_t h = 1);
    
    Tuple(uint32_t symbols, size_t h = 1);

    /**
     *  Access the \f$ i \f$-th Transposition of the Tuple.
     */
    Transposition& operator[](size_t n);

    /**
     *  @return Returns the norm (i.e. the number of transpositions).
     */
    int32_t norm() const;
    
    /**
     *  @return Returns true iff no Transposition contains the symbol 0.
     *  It is used in del2. @sa Tuple::del2.
     */
    operator bool() const;

    /** 
     *  output stream
     */
    friend std::ostream& operator<< (std::ostream& stream, const Tuple& tuple);
    
    /**
     *  Gibt an, wieviele Zykel / interne Punktierungen die Permutation \f$ \tau_h \cdot \ldots \cdot \tau_1 \cdot (1,2,\ldots,p) \f$ hat.
     */
    PermutationType permutation_type();
    
    bool monoton();                     ///< Gibt genau dann true zurueck, wenn das Tupel monoton ist.
    bool f(uint32_t i);                 ///< Die Funktion \f$ f_i \f$ fuer \f$ 1 \le i < h \f$.
    bool phi( uint32_t q, uint32_t i);   ///< Die Funtkion \f$ \Phi^q_i \f$.
    Tuple del2( uint32_t i );           ///< Der i-te horizontale Rand \f$ \partial_i^{\prime \prime} 
    
    uint32_t p;  ///< The number of symbols \f$ 1 \le p \f$ to be permuted.
private:
    std::vector< Transposition > rep;   ///< Representation eines Tuples von Transpositionen.
};

/// output stream
std::ostream& operator<< (std::ostream& stream, const Tuple& tuple);

#endif // TUPLE_HPP
