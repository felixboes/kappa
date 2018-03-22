#ifndef KAPPA_PERMUTATION_MANAGER_HPP
#define KAPPA_PERMUTATION_MANAGER_HPP

#include <algorithm>

#include "permutation.hpp"
/**
 * Class providing tools for working with permutations in form of a Permutation or a Transposition.
 * Does not contain any data.
 */
class PermutationManager {

public:

    PermutationManager();

    static bool is_trivial(const Transposition& t);

    /**
     * @return true iff the support of t1 and the support of t2 are disjoint.
     * Note that if one of t1 and t2 is trivial, then t1 and t2 have disjoint support but are not symbolwise disjoint.
     */
    static bool have_disjoint_support(const Transposition &t1, const Transposition &t2);

    /**
     * @return true iff {t1.first, t1.second} and {t2.first, t2.second} are disjoint.
     * Note that if one of t1 and t2 is trivial, then t1 and t2 have disjoint support but are not symbolwise disjoint.
     */
    static bool are_symbolwise_disjoint(const Transposition &t1, const Transposition &t2);


    /**
     * Note: if n is not in the support of t, then value_at(t,n) gives back n without checking whether n lies in the
     * domain of t as we have no domain stored for a Transposition.
     */
    static uint8_t value_at(const Transposition& t, uint8_t n);

    /**
     * Note: if n is not in the support of t, then value_at(t,n) gives back n without checking whether n lies in the
     * domain of t as we have no domain stored for a Norm2Permutation.
     */
    static uint8_t value_at(const Norm2Permutation& t, uint8_t n);

    /*
     * Note: fixed points are counted as cycles.
     */
    static uint8_t num_cycles(const Permutation& t);

    /**
     *  @brief Determines the decompositions of t1*t2 into cycles omitting fixed points.
     *
     *  Returns a vector "decomp" of vectors with:
     *  - decomp[i] is a cycle in t1*t2 written as
     *     decomp[i][0] -> decomp[i][1] ->...-> decomp[i][decomp[i].size()-1] -> decomp[i][0]
     *  - the largest symbol of decomp[i] is decomp[i][0]
     *  - the cycles are ordered descendingly in the sense that decomp[i][0]> decomp[i+1][0] for all i
     *
     *  \note easy to adapt for the product of more than two Norm2Permutations or several Transpositions
     */
    static std::vector < std::vector<uint8_t> >
    descending_cycle_decomposition_of_product(const Norm2Permutation &t1, const Norm2Permutation &t2);

    static Permutation inverse(const Permutation& t);

    /**
     * Stores the product sigma*t again in sigma.
     * Does not test whether the symbols in t are in the domain of sigma, i.e. if they are < sigma.size().
     * Note: this is way less effortful than multiplying t*sigma!
     */
    static void multiply(Permutation& sigma, const Transposition& t);

    /**
     * Stores the product sigma*t again in sigma.
     * Does not test whether the symbols in t are in the domain of sigma, i.e. if they are < sigma.size().
     * Note: this is way less effortful than computing t*sigma!
     */
    static void multiply(Permutation& sigma, const Norm2Permutation& t);

    static bool has_norm_two(const Norm2Permutation& t);

    /**
     * @brief Writes the Norm2Permutation t in norm2notation if t has norm two. Else does not change t.
     *
     * If t is a permutation of norm two, write t = t_bar * t' the factorization of t into two transpositions.
     * Write then t_bar =(i_1,i_2) with i_1 > i_2, and t'= (i_3,i_4) with
     * i_3 = height(t) > i_4. Then t = (i_1,i_2)(i_3,i_4) is the Norm2Notation of t.
     *
     * @returns true iff t has norm 2 (and thus is now in norm2notation)
     */
    static bool write_in_norm2notation(Norm2Permutation& t);

    /**
     * Note: if t is in norm2notation, then t has norm two.
     * Note: checking whether t is in norm2notation is way less effortful than calling write_in_norm2notation.
     */
    static bool is_in_norm2notation(const Norm2Permutation& t);

    /**
     * substitutes in the entries of t each 'idx' with 'replacement_idx'
     */
    static void substitute_index(Transposition &t, uint8_t idx, uint8_t replacement_idx);
    static void substitute_index(Norm2Permutation &t, uint8_t idx, uint8_t replacement_idx);

    /**
     * Removes fixed_pt in the following sense: each symbol in t larger than fixed_pt is decreased by one; each symbol
     * in t less than fixed_pt remains the same.
     * Note: fixed_pt should be a fixed point of t, i.e. should not be a symbol in t.
     */
    static void drop_fixed_point(Transposition &t, uint8_t fixed_pt);
    static void drop_fixed_point(Norm2Permutation &t, uint8_t fixed_pt);

    /**
     *  @returns the orientation signs \f$ \varepsilon_0, \ldots, \varepsilon_p \f$ of a non-degenerate cell
     *  (sigma_q,...,sigma_0) of bi-degree (p,q) with sigma = sigma_q.
     *  Note that this orientation sign indeed only depends on sigma_q.
     */
    static std::map< uint8_t, int8_t > orientation_sign_for_ehrenfried(const Permutation &sigma);

};


#endif //KAPPA_PERMUTATION_MANAGER_HPP
