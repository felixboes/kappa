#ifndef KAPPA_ALT_GRP_TUPLE_HPP
#define KAPPA_ALT_GRP_TUPLE_HPP

#include <algorithm>
#include "sym_grp_tuple.hpp"

/**
 *  A tuple of Norm2Permutations \f$ ( \tau_{num_entries} \mid \ldots \mid \tau_1 ) \f$.
 *
 *  Note: usually (eg. whenever the AltGrpTuple is a generator in the Ehrenfried complex) we are only interested in the
 *  case num_entries = h/2, where h is 2*g+m resp. 2*g+m-1 in the parallel resp. radial case. The class can be used for
 *  arbitrary num_entries anyways.
**/
class AltGrpTuple
{

public:
    friend class HashAltGrpTuple;

     /**
     *  Construct an AltGrpTuple of norm zero.
     */
    AltGrpTuple();

    /**
     *  Construct an AltGrpTuple filled with num_entries many trivial Norm2Permutations (0,0)(0,0).
     */
    explicit AltGrpTuple(size_t num_entries);

    /**
     *  Construct an AltGrpTuple on the symbols min_symbol,...,p filled with num_entries many trivial Norm2Permutations (0,0)(0,0).
     */
    AltGrpTuple(uint8_t p, size_t num_entries);

    /**
     *  Access the \f$ n \f$-th Norm2Permutation of the AltGrpTuple without bounds checked.
     */
    Norm2Permutation& operator[](size_t n);

    /**
     *  Access the \f$ n \f$-th Norm2Permutation of the AltGrpTuple with bounds checked.
     */
    Norm2Permutation& at(size_t n);

    // const version
    const Norm2Permutation& at(size_t n) const;

    /**
     *  @returns the number of Norm2Permutations in rep.
     */
    size_t num_entries() const;

    /**
     *  @returns true iff the norm is 2* num_entries, ie. iff each stored Norm2Permutation has norm two.
     */
    bool has_correct_norm() const;

    /**
     *  @return Returns true iff both AltGrpTuples are elementwise equal.
     */
    bool operator==(const AltGrpTuple& t) const;

    /**
     *  @return Returns false iff both AltGrpTuples are elementwise equal.
     */
    bool operator!=(const AltGrpTuple& t) const;

    /**
     *  Consider all AltGrpTuples as radial cells.
     */
    static void radial_case();

    /**
     *  Consider all AltGrpTuples as parallel cells.
     */
    static void parallel_case();

    static bool is_radial();

    /**
     *  @returns the minimal symbol allowed in an AltGrpTuple, i.e. 1 in the parallel case and 0 in the radial case.
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

    friend std::ostream& operator<< (std::ostream& stream, const AltGrpTuple& tuple);

    /**
     *  Computes the number of cycles of the permutation
     *   sigma_out = \f$ \tau_{num_entries} \cdot \ldots \cdot \tau_1 \cdot long_cycle(p) \f$.
     *  Fixed points are counted as cycles.
     */
    uint8_t num_cycles() const;

    /**
     * \brief Determines whether this AltGrpTuple has the correct number of cycles (including fixed points).
     * \param m number of punctures
     * \return In the parallel case: true iff the number of cycles equals m + 1;
     *         in the radial case: true iff the number of cycles equals m.
     */
    bool has_correct_num_cycles(size_t m) const;

    /**
     * Writes each Norm2Permutation \tau_i in the AltGrpTuple in Norm2Notation.
     * @returns false iff the AltGrpTuple is invalid, i.e. if some Norm2Permutation has actually not norm 2.
     */
    bool write_in_Norm2Notation();

    /**
     * @returns true iff the AltGrpTuple is fully unstable. To check this, writes the AltGrpTuple (i.e. each \tau_i) in
     * Norm2Notation.
     */
    bool fully_unstable();

    /**
     *  Applies the function \f$ f_i \f$ for \f$ 1 \le i < num_entries() \f$ and returns true iff the norm is preserved
     *  thereby.
     *  \note The parameter i must fulfill 1 <= i < num_entries()
     */
    bool f(uint8_t i);

    /**
     *  Applies the function \f$ \Phi^q_i \f$ and returns true iff the norm is preserved thereby.
     *  \note Writes the AltGrpTuple in Norm2Notation.
     *  \note The parameters must fulfill 1 <= i <= q <= num_entries()
     */
    bool phi(uint8_t q, uint8_t i);

    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} and the projection onto the fully unstable
     *  cells.
     *  Returns an empty AltGrpTuple if the boundary is degenerate or not fully unstable, the boundary AltGrpTuple in
     *  Norm2Notation otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    AltGrpTuple d_hor(uint8_t i ) const;

    /** Applies the i-th horizontal boundary  \f$ \partial_i^{\prime \prime} without the projection onto the fully
     * unstable cells.
     *  Returns an empty AltGrpTuple if the boundary is degenerate, the boundary AltGrpTuple otherwise.
     *  \note The parameter i has to fulfill 0 < i < p.
     */
    AltGrpTuple d_hor_double_complex(uint8_t i ) const;

    /**
     *  @returns the orientation signs \f$ \varepsilon_0, \ldots, \varepsilon_p \f$.
     */
    std::map< uint8_t, int8_t > orientation_sign() const;

    /**
     *  @returns the permutation sigma_out = \f$ \sigma_{num_entries} = \tau_{num_entries} \cdots \tau_1 \sigma_0 \f$,
     *  that is the permutation of the outgoing boundary curve of the (radial) slit domain.
     *
     *  \note: sigma_out_inv() is easier to compute than sigma_out().
     */
    Permutation sigma_out() const;

    /**
     *  @returns the inverse permutation of sigma_out.
     *  \note: this is easier to compute than sigma_out().
     */
    Permutation sigma_out_inv() const;


    uint8_t p; ///< The maximum of the symbols \f$ min_symbol() \le p \f$ to be permuted.
    size_t id;  ///< The index of this AltGrpTuple in the basis of the EhrComplex.

protected:

    /**
     *  @returns the cycle 0-> 1 -> 2 -> ... -> p-1 -> p -> 0.
     */
    Permutation long_cycle() const;

    Permutation long_cycle_inv() const;



    //    DATA MEMBERS

    std::vector< Norm2Permutation > rep;   ///< Representation of an AltGrpTuple by Norm2Permutations
                                           ///< \tau_{num_entries()}, ..., \tau_1, where \tau_i = rep[i-1]

    static bool radial;                    ///< true iff we consider all AltGrpTuple as radial cells.
    static uint8_t min_symbol;             ///< is set by radial_case and parallel_case accordingly.
    static uint8_t min_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.
    static uint8_t max_boundary_offset;    ///< is set by radial_case and parallel_case accordingly.

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int) ///< serialization of AltGrpTuple.
    {
        ar & p & id & rep;
    }
};


/**
 *  In order to save AltGrpTuples in a hash table (e.g. in EhrBasis) we need a function object, that hashes
 *  AltGrpTuples.
 */
class HashAltGrpTuple
{
public:
    size_t operator()( const AltGrpTuple& ) const;

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &, const unsigned int) ///< Implements the serialization.
    {
    }
};



#endif //KAPPA_ALT_GRP_TUPLE_HPP
