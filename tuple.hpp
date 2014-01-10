#ifndef TUPLE_HPP
#define TUPLE_HPP

#include <cinttypes>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include "serialization.hpp"

/**
 *  A Transposition is represented by a pair of uint8_t.
 */
typedef std::pair< uint8_t, uint8_t > Transposition;

/**
 *  A Tuple has a permutation type which stores the number of cycles of \f$ \sigma_h \f$ and the number of internal punctures.
 */ 
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
 *  We require that good Tuples use only the symbols \f$ 1 \f$ to \f$ p \f$.
**/
class Tuple
{
    friend class HashTuple;
    
public:
    typedef std::map<uint8_t, uint8_t> Permutation;  ///< Data structure to store a permutation, where a -> perm[a] if perm is a Permutation

    /**
     *  Construct a Tuple of norm h which has to be filled.
     */
    Tuple(size_t h = 1);
    
    Tuple(uint32_t symbols, size_t h);

    /**
     *  Access the \f$ i \f$-th Transposition of the Tuple.
     */
    Transposition& operator[](size_t n);

    /**
     *  @return Returns the norm (i.e. the number of transpositions).
     */
    int32_t norm() const;
    
    /**
     *  @return Returns true iff both Tuples are elementwise equal.
     */
    bool operator==(const Tuple& t) const;
    
    /**
     *  @return Returns false iff both Tuples are elementwise equal.
     */
    bool operator!=(const Tuple& t) const;
    
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
     *  Tells the number of cycles and intern punctures of the permutation \f$ \tau_h \cdot \ldots \cdot \tau_1 \cdot (1,2,\ldots,p) \f$.
     */
    PermutationType permutation_type();
    
    bool monotone();                        ///< Returns true iff the tuple is monotone.
    bool f(uint32_t i);                     ///< Applies the function \f$ f_i \f$ fuer \f$ 1 \le i < h \f$ and returns true iff the norm is preserved thereby.
    bool phi( uint32_t q, uint32_t i);      ///< Applies the function \f$ \Phi^q_i \f$ and returns true iff the norm is preserved thereby..
    Tuple d_hor( uint8_t i ) const;         ///< Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} and returns true iff the norm is preserved thereby.
    Tuple d_hor_naive( uint8_t i ) const;   ///< Different implementation of the i-th horizontal boundary.
    std::map< uint8_t, int8_t > orientation_sign( const Tuple::Permutation& sigma ) const;

    uint32_t p;  ///< The number of symbols \f$ 1 \le p \f$ to be permuted.
    uint32_t id; ///< The index of this Tuple in the basis of the MonoComplex.
private:
    Transposition& at(size_t q);                     ///< Access the q-th Transposition.
    std::vector< Transposition > rep;                ///< Representation of a tuple of transpositions.
    Permutation long_cycle() const;                  ///< Returns the cycle 1 -> 2 -> ... -> p-1 -> p -> 1.
    Permutation long_cycle_inv() const;              ///< Returns the cycle 1 -> p -> p-2 -> ... -> 2 -> 1.
    
    void print_permutation( Permutation sigma ) const;
    std::map< uint8_t, Tuple::Permutation > cycle_decomposition ( const Tuple::Permutation& sigma ) const;
    
    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization of Tuple.
    {
        ar & p & id & rep;
    }
};

/// output stream
std::ostream& operator<< (std::ostream& stream, const Tuple& tuple);

/**
 *  In order to save Tuples in a hash table (e.g. in MonoBasis) we need a function object, that hashes Tuples.
 */ 
class HashTuple
{
public:
    size_t operator()( const Tuple& ) const;
    
    friend class boost::serialization::access;
    
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) ///< Implements the serialization.
    {
    }
};

#endif // TUPLE_HPP
