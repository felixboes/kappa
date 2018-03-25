#ifndef KAPPA_PERMUTATION_HPP
#define KAPPA_PERMUTATION_HPP

#include <iostream>
#include <map>
#include <vector>

typedef std::pair< uint8_t, uint8_t > Transposition;

typedef std::pair< Transposition, Transposition > Norm2Permutation;

/**
 * @brief Class representing a permutation
 *
 * The permutation is stored in a vector of a given size, say s,
 * whereby the permutation acts on a subset of the numbers 0, ..., s-1.
 * For each symbol 0 <= i < s, its entry in the vector tells its image under the permutation.
 * Per default, each value i maps to s, marking that i does not actually belong to the permutation.
 *
 */
class Permutation
{
public:
    /**
     * Constructs an empty permutation.
     */
    Permutation();

    /**
     * Constructs a permutation of the given size and initializes each entry with the default value size.
     */
    explicit Permutation(uint8_t size);

    /**
     * Copy constructor copying the data vector from other.
     */
    Permutation(const Permutation & other);

    /**
     * Returns the data vector of this Permutation.
     */
    std::vector<uint8_t> operator()() const;

    /**
     * Returns the element the symbol i is mapped to by this Permutation without bounds checked.
     */
    uint8_t & operator[](uint8_t i);

    /**
     * Returns the element the symbol i is mapped to by this Permutation with bounds checked.
     */
    const uint8_t & at(uint8_t i) const;

    /**
     * Returns the number of elements of this permutation.
     */
    uint8_t size() const;

    /**
     *  @return Returns true iff both Permutations are elementwise equal.
     */
    bool operator==(const Permutation& t) const;

    /**
     *  @brief Determines the decompositions of the permutation into cycles including fix points.
     *  @return A map consisting of all cycles of pi with their smallest element as key
     */
    std::map< uint8_t, Permutation > cycle_decomposition () const;

    /**
     * Returns true iff the symbol i is contained in this Permutation.
     */
    bool is_contained(uint8_t i) const;

    /**
     * Returns true iff the symbol i is a fix point of this Permutation.
     */
    bool is_fix_point(uint8_t i) const;

    friend std::ostream& operator<< (std::ostream& stream, const Permutation& permutation);

    /**
     * @returns the permutation (0 -> 1 -> ... -> p -> 0).
     */
    static Permutation long_cycle(uint8_t p);

    /**
     * @returns the permutation (p -> p-1 -> ... -> 0 -> p).
     */
    static Permutation long_cycle_inv(uint8_t p);

    friend std::ostream& operator<< (std::ostream& stream, const Norm2Permutation& tau);

protected:
    std::vector<uint8_t> data; ///< stores the Permutation
    operator size_t() = delete;
};

std::ostream& operator<< (std::ostream& stream, const Norm2Permutation& tau);


#endif //KAPPA_PERMUTATION_HPP
