#ifndef KAPPA_EVENTUPLE_HPP
#define KAPPA_EVENTUPLE_HPP

#include <algorithm>
#include "sym_grp_tuple.hpp"

/**
 *  A tuple of Norm2Permutations \f$ ( \tau_{h_halves} \mid \ldots \mid \tau_1 ) \f$.
**/
class EvenTuple
{

public:
    friend class HashEvenTuple;

     /**
     *  Construct an EvenTuple of norm zero.
     */
    EvenTuple();

    /**
     *  Construct an EvenTuple filled with h_halves many trivial Norm2Permutations (0,0)(0,0).
     */
    explicit EvenTuple(size_t h_halves);

    /**
     *  Construct an EvenTuple on the symbols min_symbol,...,p filled with h_halves many trivial Norm2Permutations (0,0)(0,0).
     */
    EvenTuple(uint8_t p, size_t h_halves);

    /**
     *  Access the \f$ n \f$-th Norm2Permutation of the EvenTuple without bounds checked.
     */
    Norm2Permutation& operator[](size_t n);

    /**
     *  Access the \f$ n \f$-th Norm2Permutation of the EvenTuple with bounds checked.
     */
    Norm2Permutation& at(size_t n);

    // const version
    const Norm2Permutation& at(size_t n) const;

    /**
     *  @returns h_halves = the number of Norm2Permutations in rep.
     */
    size_t num_norm2permutations() const;

    /**
     *  @returns true iff the norm is 2* num_norm2permutation, ie. iff each stored Norm2Permutation has norm two.
     */
    bool has_correct_norm() const;

    /**
     *  @return Returns true iff both EvenTuples are elementwise equal.
     */
    bool operator==(const EvenTuple& t) const;

    /**
     *  @return Returns false iff both EvenTuples are elementwise equal.
     */
    bool operator!=(const EvenTuple& t) const;

    /**
     *  Consider all EvenTuples as radial cells.
     */
    static void radial_case();

    /**
     *  Consider all EvenTuples as parallel cells.
     */
    static void parallel_case();

    static bool is_radial();

    /**
     *  @returns the minimal symbol allowed in an EvenTuple, i.e. 1 in the parallel case and 0 in the radial case.
     */
    static uint8_t get_min_symbol();

    /**
     *  @returns the minimal face index of a (possibly non-degenerate) cell, i.e. 1 in the parallel case and 0 in the
     *  radial case.
     */
    static uint8_t get_min_boundary_offset();

    /**
     *  @returns the difference between p and the maximal face index of a (possibly non-degenerate) cell,
     *  i.e. 1 in the parallel case and 0 in the radial case.
     */
    static uint8_t get_max_boundary_offset();

    friend std::ostream& operator<< (std::ostream& stream, const EvenTuple& eventuple);

    /**
     *  Computes the number of cycles of the permutation
     *    \f$ \tau_{h_halves} \cdot \ldots \cdot \tau_1 \cdot long_cycle(p) \f$.
     *  Fixed points are counted as cycles.
     */
    uint8_t num_cycles() const;

    /**
     * \brief Determines whether this EvenTuple has the correct number of cycles (including fixed points).
     * \param m number of punctures
     * \return In the parallel case: true iff the number of cycles equals m + 1;
     *         in the radial case: true iff the number of cycles equals m.
     */
    bool has_correct_num_cycles(size_t m) const;

    /**
     * Writes each Norm2Permutation \tau_i in the EvenTuple in Norm2Notation.
     * @returns false iff the EvenTuple is invalid, i.e. if some Norm2Permutation has actually not norm 2.
     */
    bool write_in_Norm2Notation();

    /**
     * @returns true iff the EvenTuple is monotone. To check this, writes the EvenTuple (i.e. each \tau_i) in
     * Norm2Notation.
     */
    bool monotone();

    /**
     *  Applies the function \f$ f_i \f$ for \f$ 1 \le i < h_halves \f$ and returns true iff the norm is preserved
     *  thereby.
     *  \note The parameter i must fulfill 1 <= i < h_halves
     */
    bool f(uint8_t i);

    /**
     *  Applies the function \f$ \Phi^q_i \f$ and returns true iff the norm is preserved thereby.
     *  \note Does not write the EvenTuple in Norm2Notation in general.
     *  \note The parameters must fulfill 1 <= i <= q <= h_halves
     */
    bool phi(uint8_t q, uint8_t i);

    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} and the projection onto the monotonous
     *  cells.
     *  Returns an empty EvenTuple if the boundary is degenerate and monotonous, the boundary EvenTuple in Norm2Notation
     *  otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    EvenTuple d_hor(uint8_t i ) const;

    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} without the projection onto the monotonous
     *  cells.
     *  Returns an empty EvenTuple if the boundary is degenerate, the boundary EvenTuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    EvenTuple d_hor_double_complex(uint8_t i ) const;

    /**
     *  @returns the orientation signs \f$ \varepsilon_0, \ldots, \varepsilon_p \f$.
     */
    std::map< uint8_t, int8_t > orientation_sign() const;

    /**
     *  @returns the permutation \f$ \sigma_{h_halves} = \tau_{h_halves} \cdots \tau_1 \sigma_0 \f$.
     *  \note: sigma_h_halves_inv() is easier to compute than sigma_h_halves().
     */
    Permutation sigma_h_halves() const;

    /**
     *  @returns the inverse permutation of \f$ \sigma_{h_halves} \f$.
     *  \note: this is easier to compute than sigma_h_halves().
     */
    Permutation sigma_h_halves_inv() const;


    uint8_t p; ///< The maximum of the symbols \f$ min_symbol() \le p \f$ to be permuted.
    size_t id;  ///< The index of this EvenTuple in the basis of the MonoComplex.

protected:

    /**
     *  @returns the cycle 0-> 1 -> 2 -> ... -> p-1 -> p -> 0.
     */
    Permutation long_cycle() const;

    Permutation long_cycle_inv() const;



    //    DATA MEMBERS

    std::vector< Norm2Permutation > rep;   ///< Representation of an EvenTuple by Norm2Permutations
                                           ///< \tau_{h_halves}, ..., \tau_1, where \tau_i = rep[i-1]

    static bool radial;                    ///< true iff we consider all EvenTuple as radial cells.
    static uint8_t min_symbol;             ///< is set by radial_case and parallel_case accordingly.
    static uint8_t min_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.
    static uint8_t max_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int) ///< serialization of EvenTuple.
    {
        ar & p & id & rep;
    }
};


/**
 *  In order to save EvenTuples in a hash table (e.g. in MonoBasis) we need a function object, that hashes
 *  EvenTuples.
 */
class HashEvenTuple
{
public:
    size_t operator()( const EvenTuple& ) const;

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &, const unsigned int) ///< Implements the serialization.
    {
    }
};



#endif //KAPPA_EVENTUPLE_HPP
