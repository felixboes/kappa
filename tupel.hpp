#ifndef TUPEL_HPP
#define TUPEL_HPP

#include <cinttypes>
#include <map>
#include <vector>

typedef std::pair< uint8_t, uint8_t > Transposition;

/**
 *  A Tuple of Transpositions \f$ ( \tau_h \mid \ldots \mid \tau_1 ) \f$.
 *  
 *  Convention: We want an easy-to-check condition, to test if a Tuples is in a bad state.
 *  We require, that good Tuples use only the symbols \f$ 1 \f$ to \f$ p \f$.
**/
class Tupel
{
public:
    /**
     *  Construct a Tuple of norm h, which has to be filled.
     */
    Tupel(size_t h = 1);

    /**
     *  Access the \f$ i \f$-th Transposition of the Tupel.
     */
    Transposition& operator[](size_t n);

    /**
     *  @return Returns the norm (i.e. the number of transpositions).
     */
    int32_t norm() const;
    
    /**
     *  @return Returns true iff no Transposition contains the symbol 0.
     *  It is used in del2. @sa Tupel::del2.
     */
    bool valid() const;

    uint32_t p;  ///< The number of symbols \f$ 1 \le p \f$ to be permuted.
private:
    std::vector< Transposition > rep;   ///< Representation eines Tupels von Transpositionen.
};


#endif // TUPEL_HPP
